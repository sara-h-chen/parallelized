// Translate this file with
//
// g++ -O3 --std=c++11 spacebodies.c -o spacebodies
//
// Run it with
//
// ./spacebodies
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2017 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <set>
// TODO: Uncomment and test parallelization
#include <omp.h>

struct Body {
    double mass;
    double velocityX, velocityY, velocityZ;
    double positionX, positionY, positionZ;
    double forceX, forceY, forceZ;
    bool isActive;
};

const double limit = 1e-10;

double t = 0;
double tFinal = 0;
double timeStepSize = 1e-5;
bool adaptiveTimeStep = false;

int numberOfBodies = 0;
int NumInactive = 0;

Body *bodies;
std::set<std::pair<Body*, Body*>> collidedBodies;

std::ofstream videoFile;
std::ofstream bodyCountFile;


// TODO: Remove all debug statements and uncomment Paraview statements

// ----------------------------------------------------
//                  UTILITY FUNCTIONS
// ----------------------------------------------------

//void deepCopy(Body *dest, Body *src) {
//    for (int i = 0; i < numberOfBodies; ++i) {
//        dest[i].positionX = src[i].positionX;
//        dest[i].positionY = src[i].positionY;
//        dest[i].positionZ = src[i].positionZ;
//
//        dest[i].velocityX = src[i].velocityX;
//        dest[i].velocityY = src[i].velocityY;
//        dest[i].velocityZ = src[i].velocityZ;
//
//        dest[i].isActive = src[i].isActive;
//        dest[i].mass = src[i].mass;
//
//        // Set all force to 0 for calculation of next iteration
//        copyState[i].forceX = 0;
//        copyState[i].forceY = 0;
//        copyState[i].forceZ = 0;
//    }
//}

double calculateDistance(Body a, Body b) {
    return sqrt(
            (a.positionX - b.positionX) * (a.positionX - b.positionX) +
            (a.positionY - b.positionY) * (a.positionY - b.positionY) +
            (a.positionZ - b.positionZ) * (a.positionZ - b.positionZ)
    );
}

// ----------------------------------------------------
//                    SETUP FUNCTION
// ----------------------------------------------------

void setUp(int argc, char **argv) {
    numberOfBodies = (argc - 2) / 7;

    bodies = new Body[numberOfBodies];

    int readArgument = 1;

    tFinal = std::stof(argv[readArgument]);
    readArgument++;

    // Cannot be parallel because of readArgument variable
    // Process each body separately
    for (int i = 0; i < numberOfBodies; ++i) {
        // The first three numbers passed into the command line are the x variables
        // representing the position of the particle
        bodies[i].positionX = std::stof(argv[readArgument]);
        readArgument++;
        bodies[i].positionY = std::stof(argv[readArgument]);
        readArgument++;
        bodies[i].positionZ = std::stof(argv[readArgument]);
        readArgument++;

        // The next three numbers passed in represent the velocity vectors
        bodies[i].velocityX = std::stof(argv[readArgument]);
        readArgument++;
        bodies[i].velocityY = std::stof(argv[readArgument]);
        readArgument++;
        bodies[i].velocityZ = std::stof(argv[readArgument]);
        readArgument++;

        bodies[i].mass = std::stof(argv[readArgument]);
        readArgument++;
        bodies[i].isActive = true;

        if (bodies[i].mass <= 0.0) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << numberOfBodies << " bodies" << std::endl;
}

// ----------------------------------------------------
//              PARAVIEW RELATED FUNCTIONS
// ----------------------------------------------------

void openParaviewVideoFile() {
    videoFile.open("result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
              << std::endl
              << "<Collection>";
}


void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
}

/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
void printParaviewSnapshot(int counter) {
    std::stringstream filename;
    filename << "paraview/result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << numberOfBodies - NumInactive << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i = 0; i < numberOfBodies; ++i) {
        if (!bodies[i].isActive) {
            continue;
        }
        out << bodies[i].positionX
            << " "
            << bodies[i].positionY
            << " "
            << bodies[i].positionZ
            << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>"
              << std::endl;
}


// ----------------------------------------------------
//                  UPDATE FUNCTIONS
// ----------------------------------------------------

void makeForecast(Body &a, Body &b, Body &tempBodyA, Body &tempBodyB) {
    tempBodyA.positionX = a.positionX + (timeStepSize * a.velocityX);
    tempBodyB.positionX = b.positionX + (timeStepSize * b.velocityX);
    tempBodyA.positionY = a.positionY + (timeStepSize * a.velocityY);
    tempBodyB.positionY = b.positionY + (timeStepSize * b.velocityY);
    tempBodyA.positionZ = a.positionZ + (timeStepSize * a.velocityZ);
    tempBodyB.positionZ = b.positionZ + (timeStepSize * b.velocityZ);
}

void updatePosition() {
    double accelerationX, accelerationY, accelerationZ;
    double deltaT = (timeStepSize * timeStepSize) / 2;

    #pragma omp for private(accelerationX, accelerationY, accelerationZ)
    // Loop through bodies and update if active
    for (int i = 0; i < numberOfBodies; ++i) {
	// printf("Number of threads in position update: %d\n", omp_get_num_threads());
	// exit(1);
        if (bodies[i].isActive) {

            accelerationX = bodies[i].forceX / bodies[i].mass;
            accelerationY = bodies[i].forceY / bodies[i].mass;
            accelerationZ = bodies[i].forceZ / bodies[i].mass;

            bodies[i].positionX = bodies[i].positionX + (timeStepSize * bodies[i].velocityX) + (deltaT * accelerationX);
            bodies[i].positionY = bodies[i].positionY + (timeStepSize * bodies[i].velocityY) + (deltaT * accelerationY);
            bodies[i].positionZ = bodies[i].positionZ + (timeStepSize * bodies[i].velocityZ) + (deltaT * accelerationZ);

            bodies[i].velocityX = bodies[i].velocityX + (timeStepSize * accelerationX);
            bodies[i].velocityY = bodies[i].velocityY + (timeStepSize * accelerationY);
            bodies[i].velocityZ = bodies[i].velocityZ + (timeStepSize * accelerationZ);

            // Set all force to 0 for calculation of next iteration
            bodies[i].forceX = 0;
            bodies[i].forceY = 0;
            bodies[i].forceZ = 0;

            // DEBUG
//            if (bodies[i].positionX >= 0.09 && bodies[i].positionX <= 0.10) {
//                printf("\nBody %d: %7.8f  %7.8f  %7.8f", i, bodies[i].positionX, bodies[i].positionY, bodies[i].positionZ);
//            }
        }
    }
}

double scaleTimeStep(Body a, Body b, double distance) {
    // Begin breaking down
    double newTimeStep = timeStepSize;

    if (timeStepSize > limit) {
        if (distance < 0.001) {
            Body tempBodyA;
            Body tempBodyB;
            bool withinRange = true;

            makeForecast(a, b, tempBodyA, tempBodyB);
            double expectedDistance = calculateDistance(tempBodyA, tempBodyB);
//            printf("Expected distance : %7.8f \n", expectedDistance);
//            printf("Scaled distance : %7.8f \n", distance);

            // If they are getting closer together
            if (expectedDistance - distance < 0) {
                // Check if they are within range of one another
                while (withinRange) {
                    if (newTimeStep <= limit) {
                        withinRange = false;
                    } else if (expectedDistance <= limit) {
                        withinRange = false;
                    } else {
                        withinRange = (((expectedDistance - distance)/ distance) < 0.2);
                    }
                    newTimeStep = newTimeStep / 10;
                }
            }
        }
    }
    return newTimeStep;
}


void fuseBodies(Body *a, Body *b) {
    
    double combinedMass = a->mass + b->mass;
    double newVelX = ((a->mass * a->velocityX) + (b->mass * b->velocityX)) / combinedMass;
    double newVelY = ((a->mass * a->velocityY) + (b->mass * b->velocityY)) / combinedMass;
    double newVelZ = ((a->mass * a->velocityZ) + (b->mass * b->velocityZ)) / combinedMass;

    // DEBUG
    // printf("\n=====> Old body a: %5.10f, %5.10f, %5.10f, %5.10f", a->mass, a->velocityX, a->velocityY, a->velocityZ);
    // printf("\n=====> Old body b: %5.10f, %5.10f, %5.10f, %5.10f", b->mass, b->velocityX, b->velocityY, b->velocityZ);

    a->mass = combinedMass;
    a->velocityX = newVelX;
    a->velocityY = newVelY;
    a->velocityZ = newVelZ;

    // DEBUG
    // printf("\n=====> New combined body : %5.10f, %5.10f, %5.10f, %5.10f", a->mass, a->velocityX, a->velocityY, a->velocityZ);

    b->isActive = false;
    // Protect this variable against concurrent access
    #pragma omp atomic
    NumInactive += 1;
}

void addForce(Body *a, Body *b, double distance) {
    double massDistance = a->mass * b->mass / (distance * distance * distance);
    // Calculate the force between body and the others
    a->forceX += (b->positionX - a->positionX) * massDistance;
    a->forceY += (b->positionY - a->positionY) * massDistance;
    a->forceZ += (b->positionZ - a->positionZ) * massDistance;
}

// Calculate the distance; add force to the particle based on this distance
void calculateEffect(int a_index, int b_index) {
    Body *a = &bodies[a_index];
    Body *b = &bodies[b_index];

    // Current state
    double distance = calculateDistance(*a, *b);
    // DEBUG
//    if (distance < 1e-5) {
//        printf("Distance: %7.64f \n", distance);
//        printf("\nBody %d: %7.64f  %7.64f  %7.64f", a_index, bodies[a_index].positionX, bodies[a_index].positionY, bodies[a_index].positionZ);
//        printf("\nBody %d: %7.64f  %7.64f  %7.64f", b_index, bodies[b_index].positionX, bodies[b_index].positionY, bodies[b_index].positionZ);
//    }


    // Apply smallest time step to global
    if (adaptiveTimeStep) {
        double thisTimeStep = scaleTimeStep(*a, *b, distance);
        if (thisTimeStep < timeStepSize) {
            #pragma omp critical
            timeStepSize = thisTimeStep;
        }
    }

//    if (adaptiveTimeStep) {
//        timeStepSize = scaleTimeStep(*a, *b, distance);
        // DEBUG
//        if (distance <= 1e-8) {
//            printf("adaptive time step: %7.30f ", timeStepSize);
//        }
//    }

    if (distance > 1e-8) {
        addForce(a, b, distance);
        addForce(b, a, distance);
	return;
    } 
        // DEBUG
//        printf("\n -------------- Bodies %d and %d should collide with distance : %5.40f\n", a_index, b_index, distance);
//        printf("\nBody %d: %7.64f  %7.64f  %7.64f", a_index, bodies[a_index].positionX, bodies[a_index].positionY, bodies[a_index].positionZ);
//        printf("\nBody %d: %7.64f  %7.64f  %7.64f", b_index, bodies[b_index].positionX, bodies[b_index].positionY, bodies[b_index].positionZ);
}

// Part 3: Make the time step change according to how close the bodies are to one another so that the particles don't just pass through each other
void updateBodies() {
    // Step 1.1: All bodies interact and move
    #pragma omp for 
    for (int i = 0; i < numberOfBodies; ++i) {
        if (bodies[i].isActive) {
            for (int j = i + 1; j < numberOfBodies; ++j) {
                if (bodies[j].isActive) {
                    calculateEffect(i, j);
                }
            }
        }
    }

    updatePosition();

    // Increase time
    // Guarantee thread safety
    #pragma omp atomic
    t += timeStepSize;
}

void addToSet(Body *a, Body *b) {
    std::pair<Body*, Body*> first_perm(*a, *b);
    std::pair<Body*, Body*> second_perm(*b, *a);
    collidedBodies.insert(first_perm, second_perm); 
}


// ----------------------------------------------------
//                RANDOM BODY GENERATOR
// ----------------------------------------------------

void createRandomBodies(int noOfBodies) {
    std::random_device rd;
    std::default_random_engine e2(rd());
    std::uniform_real_distribution<> pos_dist(-1, 1);
    std::uniform_real_distribution<> vel_dist(-10, 10);
    std::uniform_real_distribution<> mass_dist(0, 1);

    bodies = new Body[noOfBodies];

    #pragma omp parallel for
    for (int i = 0; i < noOfBodies; ++i) {   
        bodies[i].positionX = pos_dist(e2);
        bodies[i].positionY = pos_dist(e2);
        bodies[i].positionZ = pos_dist(e2);

        bodies[i].velocityX = vel_dist(e2);
        bodies[i].velocityY = vel_dist(e2);
        bodies[i].velocityZ = vel_dist(e2);

        bodies[i].mass = mass_dist(e2);
        bodies[i].isActive = true;
    }

    // DEBUG
//    std::cout << "created random setup with " << noOfBodies << " bodies" << std::endl;
}


// ----------------------------------------------------
//                 COMMAND LINE PARSER
// ----------------------------------------------------
//  If the -r flag is used, create random bodies.
//  Must be followed by the number of bodies.
// ----------------------------------------------------

bool checkFlag(char** begin, char** end, const std::string &option)
{
    return std::find(begin, end, option) != end;
}

char* getCmdOption(char** begin, char** end, const std::string &option)
{
    char** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

// ----------------------------------------------------
//                     MAIN METHOD
// ----------------------------------------------------

int main(int argc, char **argv) {

    int nProcessors = omp_get_max_threads();
    omp_set_num_threads(nProcessors);

    clock_t tStart;
    tStart = clock();

    // TODO: Remove later when not required
    // Run smaller parallel simulation

    // Check if create bodies
    if(checkFlag(argv, argv + argc, "-r")) {
        // Set time step
        tFinal = 0.1;

        // Get number of bodies from command line
        char* cmdOption = getCmdOption(argv, argv + argc, "-r");
        numberOfBodies = atoi(cmdOption);
        createRandomBodies(numberOfBodies);
        
	cmdOption = getCmdOption(argv, argv+argc, "-p");
        std::stringstream filename;
        filename << "./body_count_ten_thousand_" << cmdOption << ".csv";
	bodyCountFile.open(filename.str().c_str());

    } else {
        // Insufficient args
        if (argc == 1) {
            std::cerr << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m"
                      << std::endl
                      << std::endl
                      << "Examples:" << std::endl
                      << "100.0   0 0 0   1.0 0 0   1.0 \t One body moving form the coordinate system's centre along x axis with speed 1"
                      << std::endl
                      << "100.0   0 0 0   1.0 0 0   1.0     0 1.0 0   1.0 0 0   1.0 \t One spiralling around the other one"
                      << std::endl
                      << "100.0   3.0 0 0   0 1.0 0   0.4     0 0 0   0 0 0   0.2     2.0 0 0   0 0 0   1.0 \t Three body setup from first lecture"
                      << std::endl
                      << std::endl;
            return -1;
        }
            // Mismatched args
        else if ((argc - 2) % 7 != 0) {
            std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)"
                      << std::endl;
            return -2;
        }

        setUp(argc, argv);
    }

    // openParaviewVideoFile();
    // printParaviewSnapshot(0);

    int currentTimeSteps = 0;
    const int plotEveryKthStep = 100;
    // OpenMP keeps pools of threads alive for a bit so if generating random bodies;
    // new threads may not be spawned
    #pragma omp parallel num_threads(4)
    while (t <= tFinal) {
            // if (currentTimeSteps % plotEveryKthStep == 0) {

            // DEBUG
//             std::cout << "Going into snapshot " << currentTimeSteps/plotEveryKthStep << std::endl;

            // Print number of bodies
            bodyCountFile << currentTimeSteps << ", " << std::to_string(numberOfBodies - NumInactive) << std::endl;

            // printParaviewSnapshot(currentTimeSteps/plotEveryKthStep);
        // }
        updateBodies();

        // Continue with only one thread
        // Ensures that a particle is not fused by a thread and another thread tries with the old copy of the particle
        fuseBodies(a, b);
	#pragma omp atomic
        currentTimeSteps++;
    }

    // closeParaviewVideoFile();

    tStart = clock() - tStart;
    double time_taken = ((double)tStart)/CLOCKS_PER_SEC;
    printf("Time taken: %f \n", time_taken);
//    std::cout << time_taken << std::endl;

    bodyCountFile.close();
    return 0;
}
