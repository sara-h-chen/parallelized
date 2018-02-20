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

// NOTE: Please ignore any commented out DEBUG sections

#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <map>
#include <omp.h>

struct Body {
    double mass;
    double velocityX, velocityY, velocityZ;
    double positionX, positionY, positionZ;
    double forceX, forceY, forceZ;
    bool isActive;
};

// PART 1 CALCULATION
// int shrinkTime = 64;
// const double defaultTime = 1e-5 / shrinkTime;

const double limit = 1e-10;
const double defaultTime = 1e-5;
const int seed = 1337;

double t = 0;
double tFinal = 0;
double timeStepSize = 1e-5;
bool adaptiveTimeStep = false;

int numberOfBodies = 0;
int NumInactive = 0;

Body *bodies;

// For scaling the adaptive time step
std::pair<Body, Body> closestPair;
double closestDistance = std::numeric_limits<double>::max();

// Prepare data structure to note collision
std::map<Body*, Body*> collidedBodies;

std::ofstream videoFile;
std::ofstream bodyCountFile;


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

// SAME = 0; NOT SAME = 1 
// TRUE = 1; FALSE = 0
void addToMap(Body *a, Body *b) {
    // Only add if both bodies are not present in the map
    #pragma omp critical (Map)
    if (((collidedBodies.find(a)->second) - b) && ((collidedBodies.find(b)->second) - a)) {
	// DEBUG
	// printf("A: %p\n",  (void *)collidedBodies.find(a)->second);
	// printf("B: %p\n",  (void *)collidedBodies.find(b)->second);
	// std::cout<< (collidedBodies.find(a)->second) - b << std::endl;
        collidedBodies.insert(std::make_pair(a,b));
    }
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

	bodies[i].forceX = 0.0;
	bodies[i].forceY = 0.0;
	bodies[i].forceZ = 0.0;

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
    videoFile.open("paraview/result.pvd");
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
    // Count number of Active Bodies
    // int doCount = 0;
    // for (int i = 0; i < numberOfBodies; ++i) {
    //     if(bodies[i].isActive) {
    //         doCount += 1;
    //  }
    // }

    std::stringstream filename;
    filename << "paraview/result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << (numberOfBodies - NumInactive) << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    out << std::fixed << std::setprecision(5); 
    for (int i = 0; i < numberOfBodies; ++i) {
        if (bodies[i].isActive) {
            out << bodies[i].positionX
                << " "
                << bodies[i].positionY
                << " "
                << bodies[i].positionZ
                << " \n";
        }
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

    #pragma omp parallel for private(accelerationX, accelerationY, accelerationZ)
    // Loop through bodies and update if active
    for (int i = 0; i < numberOfBodies; ++i) {
	// DEBUG
	// printf("Number of bodies: %d\n", numberOfBodies);
	// printf("Number of threads in position update: %d\n", omp_get_num_threads());
        // printf("Number of threads in position update: %d\n", omp_get_thread_num());
	// exit(1);
        if (bodies[i].isActive) {

            accelerationX = bodies[i].forceX / bodies[i].mass;
            accelerationY = bodies[i].forceY / bodies[i].mass;
            accelerationZ = bodies[i].forceZ / bodies[i].mass;

            bodies[i].positionX = bodies[i].positionX + (timeStepSize * bodies[i].velocityX) + (deltaT * accelerationX);
            bodies[i].positionY = bodies[i].positionY + (timeStepSize * bodies[i].velocityY) + (deltaT * accelerationY);
	    // DEBUG
	    // printf("y=%.30f, v=%f, dta=%f  \n", bodies[i].positionY, bodies[i].velocityY, deltaT * accelerationY);
            bodies[i].positionZ = bodies[i].positionZ + (timeStepSize * bodies[i].velocityZ) + (deltaT * accelerationZ);

            bodies[i].velocityX = bodies[i].velocityX + (timeStepSize * accelerationX);
            bodies[i].velocityY = bodies[i].velocityY + (timeStepSize * accelerationY);
            bodies[i].velocityZ = bodies[i].velocityZ + (timeStepSize * accelerationZ);

            // Set all force to 0 for calculation of next iteration
            bodies[i].forceX = 0;
            bodies[i].forceY = 0;
            bodies[i].forceZ = 0;

            // DEBUG
            // if (bodies[i].positionX >= 0.09 && bodies[i].positionX <= 0.10) {
            //     printf("\nBody %d: %7.8f  %7.8f  %7.8f", i, bodies[i].positionX, bodies[i].positionY, bodies[i].positionZ);
            // }
        }
    }
}

double scaleTimeStep(Body a, Body b, double distance) {
    // Begin scaling down from initial value
    double newTimeStep = defaultTime;

    if (timeStepSize > limit) {
        if (distance < 0.001) {
            Body tempBodyA;
            Body tempBodyB;
            bool withinRange = true;

            makeForecast(a, b, tempBodyA, tempBodyB);
            double expectedDistance = calculateDistance(tempBodyA, tempBodyB);
	    // DEBUG
            // printf("Expected distance : %7.8f \n", expectedDistance);
            // printf("Scaled distance : %7.8f \n", distance);

            // If they are getting closer together
            if (expectedDistance - distance < 0) {
                // Check if they are within range of one another
                while (withinRange) {
                    if (newTimeStep <= 1e-11) {
                        withinRange = false;
		    // Check if the distance is too small
                    } else if (expectedDistance <= limit) {
                        withinRange = false;
		    // Check that they are moving closer together and not going too far
                    } else {
                        withinRange = (((expectedDistance - distance)/ distance) < 0.2);
                    }
                    newTimeStep = newTimeStep / 2;
                }
            }
        }
    }
    // DEBUG
    // printf("Time step: %3.15f\n", newTimeStep);
    return newTimeStep;
}


void fuseBodies(Body *a, Body *b) {
    if (a->isActive && b->isActive) {
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
        NumInactive += 1;
    }
}

void addForce(Body *a, Body *b, double distance) {
    double massDistance = (a->mass * b->mass) / (distance * distance * distance);
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
    // if (distance < 1e-5) {
    //     printf("Small distance: %7.64f \n", distance);
    //     printf("\nBody %d: %7.64f  %7.64f  %7.64f", a_index, bodies[a_index].positionX, bodies[a_index].positionY, bodies[a_index].positionZ);
    //     printf("\nBody %d: %7.64f  %7.64f  %7.64f", b_index, bodies[b_index].positionX, bodies[b_index].positionY, bodies[b_index].positionZ);
    // }

    if (adaptiveTimeStep) {
	if (distance < closestDistance) {
	    closestDistance = distance;
	    closestPair = std::make_pair(*a, *b);
	}
    }

    if (distance > 1e-8) {
        addForce(a, b, distance);
        addForce(b, a, distance);
    } else {
        // DEBUG
        // printf("\n -------------- Bodies %d and %d should collide with distance : %5.40f\n", a_index, b_index, distance);
        addToMap(a, b);
        // DEBUG
        // printf("\nBody %d: %7.64f  %7.64f  %7.64f", a_index, bodies[a_index].positionX, bodies[a_index].positionY, bodies[a_index].positionZ);
        // printf("\nBody %d: %7.64f  %7.64f  %7.64f", b_index, bodies[b_index].positionX, bodies[b_index].positionY, bodies[b_index].positionZ);
	// exit(-2);
    } 
}

// Part 3: Make the time step change according to how close the bodies are to one another so that the particles don't just pass through each other
void updateBodies() {
    
    // Clock statements here are used to time the parallel
    // and serial parts of the code
    clock_t timeStart;
    timeStart = clock();

    // Step 1.1: All bodies interact and move
    #pragma omp parallel 
    #pragma omp for nowait
    for (int i = 0; i < numberOfBodies; ++i) {
        if (bodies[i].isActive) {
            for (int j = i + 1; j < numberOfBodies; ++j) {
                if (bodies[j].isActive) {
                    calculateEffect(i, j);
                }
            }
        }
    }

    timeStart = clock() - timeStart;
    double time_count = ((double)timeStart)/CLOCKS_PER_SEC;
    // printf("Time taken for parallel force calculation: %5.15f", time_count);
    timeStart = clock();
    
    // Continue with only one thread
    if (adaptiveTimeStep) {	
        timeStepSize = scaleTimeStep(closestPair.first, closestPair.second, closestDistance);
	// DEBUG
	// printf("Scaled time step: %.15f", timeStepSize);
    }
    
    updatePosition();

    timeStart = clock() - timeStart;
    time_count = ((double)timeStart)/CLOCKS_PER_SEC;
    // printf("Time taken for parallel update: %5.15f", time_count);
    timeStart = clock();

    // Fuse the bodies serially to ensure that the process is threadsafe
    for (auto const &x: collidedBodies) { 
        fuseBodies(x.first, x.second);
    }
    collidedBodies.clear();
        
    // Increase time
    t += timeStepSize;
    timeStepSize = defaultTime;
    timeStart = clock() - timeStart;
    time_count = ((double)timeStart)/CLOCKS_PER_SEC;
}

// ----------------------------------------------------
//                RANDOM BODY GENERATOR
// ----------------------------------------------------

void createRandomBodies(int noOfBodies) {
    std::random_device rd;
    std::default_random_engine e2(rd());
    // Seeded random number generator to test similarity between
    // serial and parallel code
    // std::default_random_engine e2(seed);
    
    // Generate random numbers within this range
    std::uniform_real_distribution<> pos_dist(-1e-3, 1e-3);
    std::uniform_real_distribution<> vel_dist(-9, 9);
    std::uniform_real_distribution<> mass_dist(1e-9, 1e-10);

    bodies = new Body[noOfBodies];

    for (int i = 0; i < noOfBodies; ++i) {   
        bodies[i].positionX = pos_dist(e2);
        bodies[i].positionY = pos_dist(e2);
        bodies[i].positionZ = pos_dist(e2);

        bodies[i].velocityX = 0.0; 
        bodies[i].velocityY = 0.0; 
        bodies[i].velocityZ = 0.0;

        bodies[i].mass = mass_dist(e2);
        bodies[i].isActive = true;

	bodies[i].forceX = 0.0;
	bodies[i].forceY = 0.0;
	bodies[i].forceZ = 0.0;

	// DEBUG
	// printf("Body %d: %3.5f x, %3.5f y, %3.5f z", i, bodies[i].positionX, bodies[i].positionY, bodies[i].positionZ);
    }
    // DEBUG
    // std::cout << "created random setup with " << noOfBodies << " bodies" << std::endl;
}


// ----------------------------------------------------
//                 COMMAND LINE PARSER
// ----------------------------------------------------
//  If the -r flag is used, create random bodies.
//  Must be followed by the number of bodies.
//  The -p flag is to add a suffix to filenames.
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

    // Set number of threads explicitly or maximize
    // int nProcessors = omp_get_max_threads();
    omp_set_num_threads(4);

    clock_t tStart;
    tStart = clock();

    // Check if to generate random bodies
    if(checkFlag(argv, argv + argc, "-r")) {
        // Set length of simulation 
        tFinal = 0.1;

        // Get number of bodies from command line
        char* cmdOption = getCmdOption(argv, argv + argc, "-r");
        numberOfBodies = atoi(cmdOption);
        createRandomBodies(numberOfBodies);
        
	// Get suffix for file to output checkpoints containing
	// time steps taken, bodies, clock time
	cmdOption = getCmdOption(argv, argv+argc, "-p");
        std::stringstream filename;
        filename << "./amdahl_" << cmdOption << ".csv";
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

	// DEBUG
	bodyCountFile.open("calc.csv");
        setUp(argc, argv);
    }

    // UNCOMMENT TO PRINT PARAVIEW FILES
    // openParaviewVideoFile();
    // printParaviewSnapshot(0);

    int currentTimeSteps = 0;
    const int plotEveryKthStep = 100;
    while (t <= tFinal) {
         if (currentTimeSteps % plotEveryKthStep == 0) {

            // DEBUG
            // std::cout << "Going into snapshot " << currentTimeSteps/plotEveryKthStep << std::endl;

	    clock_t tTime;
	    tTime = clock();
	    tTime = clock() - tStart;
	    double time_count = ((double)tTime)/CLOCKS_PER_SEC;
            // Print time steps, clock time, no. of bodies to file ending with -p arg
            bodyCountFile << "time units: " << t << ", timesteps: " << currentTimeSteps << " " << std::to_string(time_count) << ", bodies left: " << std::to_string(numberOfBodies - NumInactive) << std::endl;

	    // UNCOMMENT TO PRINT PARAVIEW FILES
            // printParaviewSnapshot(currentTimeSteps/plotEveryKthStep);
        }
        updateBodies();
	currentTimeSteps++;
    }

    // UNCOMMENT TO PRINT PARAVIEW FILES
    // closeParaviewVideoFile();

    // Print time taken to run entire simulation
    tStart = clock() - tStart;
    double time_taken = ((double)tStart)/CLOCKS_PER_SEC;
    printf("Time taken: %f \n", time_taken);

    bodyCountFile.close();
    return 0;
}
