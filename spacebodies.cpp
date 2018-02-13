// Translate this file with
//
// g++ -O3 --std=c++11 spacebodies.cpp -o spacebodies
//
// Run it with
//
// ./spaceboddies
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2017 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <ctime>

#include <string>
#include <math.h>
#include <unordered_map>
#include <fstream>
#include <omp.h>
#include <time.h>
#include <iomanip>

// CONFIGS
bool adaptiveTimeStepCheck = false; // If true, then use adaptive timestep
bool useParaview = true; // if true, write paraview related files.
bool runParallel = true; // if set to parallel; then we run it in series.

double defaultTimeStepSize = 0.01;
double smallSizeLimit = 1e-8; // If variables are smaller than this then it might as well be zero.

int numberOfIterations = 20; // When using the collision iteration (for different timesteps)
bool isCsvCollisionWrite = false; // if set to true, it will generate a csv of two bodies and data to show their collision.
bool collisionIterate = false; // iterates through multiple rounds, halving the timestep size as it goes.

bool isCsvBodyCountWrite = true; // if true, will write a csv that counts the number of bodies over time.

// writer variables
bool printBodiesInfo = false; // prints the bodies and their data
bool printTimestampInfo = true;// print timestamp if true

double referenceDistance = -0.45; // If measuring error, use this value as the reference!

// Tools to manage random spacebodies
// (This is utilised by initiating spacebodies without any parameters.)
bool RandomBodies = true;
int NumberOfRandomBodies = 10;
double randomSimTFinal = 10.0;
double fMin = -1.00; // random double min
double fMax = 1.00; // random double max

double seed = 12321321;

// SETUP VARIABLES-------------------------------

  bool isCollided = false;
  double t = 0;
  double tFinal = 0;
  int NumberOfBodies = 0;
  double currentTimestep = 0;
  double numberOfCores = 1;
  int timeStepsSinceLastPlot = 0;
  std::ofstream videoFile;
  clock_t tStart = clock();
  int initialNumberOfBodies = 0;

// CSV WRITER HANDLERS --------------------------

  std::ofstream csvBodyCountFile ("bodycount.csv");
  std::ofstream csvCollisionFile ("collision.csv");

  // helper functions to write to csv.
  void countWrite(std::string text){
    if (isCsvBodyCountWrite){
      csvBodyCountFile << text;
    }
  }
  void collisionWrite(std::string text){
    if (isCsvCollisionWrite){
      csvCollisionFile << text;
    }
  }

// BODY CLASS -----------------------------------

  class Body {
    public:
      double force[3]; //body force
      double x[3];   // displacement of body
      double v[3];   // velocity of body
      double mass;    // body mass

      void reset(){
        for (int i = 0; i < 3; i++){
          force[i] = 0;
          x[i] = 0;
          v[i] = 0;
        }
        mass = 0;
      }

      void capForceLimit(){
        // If the bodies force is small enough then we set it to zero.
        for (int k=0; k<3; k++){if (force[k] < 1e-10){force[k] = 0;}}
      }

      void printForce(){
        printf("Force: %+020.15f, %+020.15f, %+020.15f \n", force[0],force[1],force[2]);
      }

      void resetForce(){
        for (int i = 0; i < 3; i++) {force[i] = 0;}
      }

      void print(int bodyID){
        // literally prints body statistics.
        printf("Body %4d: %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f \n", bodyID, x[0], x[1], x[2], v[0], v[1], v[2], mass);
      }
  };

// PARAVIEW FUNCTIONS ---------------------------

  /*
  * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
  */

  void printParaviewSnapshot(int counter, Body * bodies) {
    std::stringstream filename;
    filename << "paraview/result-" << counter <<  ".vtp";
    std::ofstream out( filename.str().c_str() );
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i=0; i<NumberOfBodies; i++) {
      out << bodies[i].x[0]
          << " "
          << bodies[i].x[1]
          << " "
          << bodies[i].x[2]
          << " \n";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>"  << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
  }

  void openParaviewVideoFile() {
    videoFile.open( "paraview/result.pvd" );
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "<Collection>";
  }

  void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
  }

// COLLISION HANDLERS ---------------------------

  void collisionDebug(Body a, Body b){
    if (isCsvCollisionWrite){
      double pos = a.x[0]-referenceDistance;
      
      std::ostringstream strs;
      strs << pos;
      std::string errStr = strs.str();

      std::ostringstream str2;
      str2 << currentTimestep;
      std::string ts = str2.str();

      std::ostringstream str3;
      str3 << currentTimestep / pos;
      std::string t3 = str3.str();

      // std::cout << "COLLISION @ t="<< t<<", Timestep: "<< currentTimestep << ": Error: " << pos << ", Location: " << a.x[0] << "; " << b.x[0] << ", Ratio: " << (currentTimestep/pos) <<", Steps: "<< timeStepsSinceLastPlot << ", Range: " << abs(b.x[0] - a.x[0]) << std::endl;

      std::cout <<currentTimestep<<"," << pos << "," << a.x[0] << "," << b.x[0] << "," << (currentTimestep/pos) <<","<< timeStepsSinceLastPlot << "," << abs(b.x[0] - a.x[0]) << std::endl;
      
      collisionWrite(ts+ ", " + errStr + ","+t3+"\n");
    }
  }


  // Fuses two bodies together.
  Body joinBodies(Body x, Body y){
    Body z;
    // combine mass.
    z.mass = x.mass + y.mass;
    for (int i; i < 3; i++){
      // get position of one of the colluded bodies
      z.x[i] = (x.x[i] + y.x[i]) / 2;
      // calculate velocity of new item. (Mass averaged velocities)
      z.v[i] = (x.v[i]/x.mass) + (y.v[i]/y.mass);
    }
    return z;
  }

  // Checks 
  void checkCollision(Body b[], int positions[][2], int collisionPairCount){
    for (int i = 0; i < collisionPairCount; i++){
       collisionDebug(b[positions[i][0]], b[positions[i][1]]);
    }
    std::unordered_map<std::string,bool> used_bodies;

    // initiate new pos of collided bodies.
    int replacementBodies[NumberOfBodies];
    for (int i = 0; i < NumberOfBodies; i++){
      replacementBodies[i] = -1;
    }

    while (collisionPairCount > 0){
      // out of the two items, which one is the smallest one? (doesn't matter rly)
      int u = positions[collisionPairCount-1][0];
      int v = positions[collisionPairCount-1][1];

      // if u has been used prior, replace u with its offset value (cus its been)
      // made already.
      if (replacementBodies[u] != -1){u = replacementBodies[u];}
      if (replacementBodies[v] != -1){v = replacementBodies[v];}

      // disgusting swap of u and v.
      int k = 0;
      if (u > v){k = u; u = v; v = k;}

      // if they're the same then do nothing..
      if (u - v != 0){
        auto identifier = std::to_string(u) + "_" + std::to_string(v);

        if (used_bodies[identifier] != true){
          // replace body at smaller index with the colluded bodies.
          b[u] = joinBodies(b[u], b[v]);
          // swap higher index with the last item in the list
          b[v] = b[NumberOfBodies-1];
          // add the replacement for v to go to u.
          replacementBodies[v] = u;
          used_bodies[identifier] = true;
          // since we've merged bodies, reduce the body count.
          NumberOfBodies--;
        }
      }
      // update iterations.
      collisionPairCount--;
    }

  }

// CORE FUNCTIONS -------------------------------

  // read arguments
  void setUp(int argc, char** argv, Body* b) {
    int readArgument = 1;
    // Interprets a floating point value in a string str
    tFinal = std::stof(argv[readArgument]); readArgument++;
    for (int i=0; i<NumberOfBodies; i++) {
      // get position values
      b[i].x[0] = std::stof(argv[readArgument]); readArgument++;
      b[i].x[1] = std::stof(argv[readArgument]); readArgument++;
      b[i].x[2] = std::stof(argv[readArgument]); readArgument++;
      // get velocity values
      b[i].v[0] = std::stof(argv[readArgument]); readArgument++;
      b[i].v[1] = std::stof(argv[readArgument]); readArgument++;
      b[i].v[2] = std::stof(argv[readArgument]); readArgument++;
      // get mass value
      b[i].mass = std::stof(argv[readArgument]); readArgument++;

      if (b[i].mass <= 0.0) {
        std::cout << "invalid mass for body " << i << std::endl;
        exit(-2);
      }
      // printf("Body %5.0d \t %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n", i, b[i].x[0], b[i].x[1], b[i].x[2], b[i].v[0], b[i].v[1], b[i].v[2], b[i].mass);
    }
    // std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  }

  // manipulate the timestep (only called when adaptive is enabled)
  double manipulateTimestep(double beforeTS, Body a, Body b, double currentDistance){
    double timestep = beforeTS;
    // if its too small don't change the timestep or we'll get nowhere.
    if (timestep > smallSizeLimit){
      // if the smallest distance is small enough we'll update the timestep to accommodate it.
      if (currentDistance < 0.01){
        // placeholder values to estimate future positions
        double newAX[3];
        double newBX[3];
        double newDist;
        bool areTheyGettingCloser = true;
        double microStep = 1e-7;

        // lets see if the bodies are getting closer or going away.
        for (int i = 0; i < 3; i++){
          newAX[i] = (microStep * a.v[i]) + a.x[i];
          newBX[i] = (microStep * b.v[i]) + b.x[i];
        }
        newDist = sqrt(
          (newBX[0]-newAX[0]) * (newBX[0]-newAX[0]) +
          (newBX[1]-newAX[1]) * (newBX[1]-newAX[1]) +
          (newBX[2]-newAX[2]) * (newBX[2]-newAX[2])
        );
        if((newDist - currentDistance) * 1e7 > 0){
          // if they're going further away then we don't need to bother about it.
          areTheyGettingCloser = false;
        }

        // comprehensible code
        bool badDistanceRatio = areTheyGettingCloser;
  
        // if they are within the region, then we adjust the timestep
        // so that the future step would be a fraction of the distance
        // they just covered.
        while (badDistanceRatio){
          // reduce the timestep.
          timestep = timestep / 10;
          // if this timestep is small enough; then breakout
          if (timestep < smallSizeLimit){
            break;
          }
          // calculate the next distance.
          for (int i = 0; i < 3; i++){
            newAX[i] = (timestep * a.v[i]) + a.x[i];
            newBX[i] = (timestep * b.v[i]) + b.x[i];
          }
          newDist = sqrt(
            (newBX[0]-newAX[0]) * (newBX[0]-newAX[0]) +
            (newBX[1]-newAX[1]) * (newBX[1]-newAX[1]) +
            (newBX[2]-newAX[2]) * (newBX[2]-newAX[2])
          );
          // they could actually collide at this stage; lets check for that first.
          if (newDist <= smallSizeLimit){
            // break the loop, we have a perfect score.
            badDistanceRatio = false;
          } else {
            // we want the new distance to be around 1/3 to 1/2 of the distance they will cover.
            badDistanceRatio = (((newDist-currentDistance) / currentDistance) < 0.2);
          }
        }
      }

      // vanilla adaptive timestepping method; for use when everything else fails
      // double EPS = 1e-8;
      // for (int i = 0; i < 3; i++){

      //   if (a.v[i] > 0 and (abs(timestep * a.force[0]/a.mass)/a.v[0] > EPS)){
      //     timestep = timestep/2;
      //   } 
      // }
    }
    return timestep;
  }

  // function that updates the positions of the particles in space
  void updateBodies(Body* bodies) {
    int i,j,k; //incrementer variables
    currentTimestep = defaultTimeStepSize;
    if (printBodiesInfo){printf("\n\n");}
    if (printTimestampInfo){printf ("T: %012.10f, # Bodies: %1.0d \n", t, NumberOfBodies);}

    // initiate positions of the shortest body positions
    double closestDistance = 999999;
    int closestPair1 = 0;
    int closestPair2 = 0;
    // initiate values we'll use to determine collided bodies
    int collisions[NumberOfBodies][2];
    int collisionPairCount = 0;
    int BeforeBodyCount = NumberOfBodies;

    #pragma omp parallel default(none) private(j,k) shared(i, isCollided, smallSizeLimit, adaptiveTimeStepCheck, BeforeBodyCount, printBodiesInfo, t, bodies, currentTimestep, closestDistance, closestPair1, closestPair2, collisions, collisionPairCount, NumberOfBodies)
    // #pragma omp parallel
    {
      
      // Calculate distances and forces for each body
      // Note that we can run this in SIMD since none of the 
      // iterations depends on the previous ones
      #pragma omp for schedule(guided)
      for (i=0; i<NumberOfBodies; i++) {
        // print bodies if enabled
        if (printBodiesInfo){bodies[i].print(i);}
        // only calculate half the calculations!
        for (j=0; j<i; j++) {
          if (i != j){ // make sure it doesn't interact with itself.
            // calculate distance
            
            double distance = sqrt(
              (bodies[i].x[0]-bodies[j].x[0]) * (bodies[i].x[0]-bodies[j].x[0]) +
              (bodies[i].x[1]-bodies[j].x[1]) * (bodies[i].x[1]-bodies[j].x[1]) +
              (bodies[i].x[2]-bodies[j].x[2]) * (bodies[i].x[2]-bodies[j].x[2])
            );
            // have they collided?
            if (distance <= smallSizeLimit){
              isCollided = true;
              collisions[collisionPairCount][0] = i;
              collisions[collisionPairCount][1] = j;
              collisionPairCount++;
            } else {
              // check if this is the closest body to date
              if (closestDistance > distance){
                closestDistance = distance;
                closestPair1 = i;
                closestPair2 = j;
              }
            }
            // calculate combined mass
            double calc4 = bodies[i].mass*bodies[j].mass/distance/distance/distance;

            // update force for both i and j (to reduce the number of ops)
            for (k = 0; k < 3; k++){
              bodies[i].force[k] += (bodies[j].x[k]-bodies[i].x[k]) * calc4;
              bodies[j].force[k] += (bodies[i].x[k]-bodies[j].x[k]) * calc4;
            }
          }
        }
      }
    // }
      #pragma omp barrier

      #pragma omp single
      {
        // There isn't really a point on parallelising this
        if (adaptiveTimeStepCheck == true){
          // update timestep if needed to accommodate the smallest distance
          currentTimestep = manipulateTimestep(currentTimestep, bodies[closestPair2], bodies[closestPair1], closestDistance);
        }

        // initialising multiple threads will probably take longer.
        for (i=0; i<NumberOfBodies; i++){
          // update position and velocity of body
          for (j=0; j<3; j++){
            // update the position and velocity of the body.
            bodies[i].x[j] = bodies[i].x[j] + (currentTimestep * bodies[i].v[j]);
            bodies[i].v[j] = bodies[i].v[j] + (currentTimestep * (bodies[i].force[j] / bodies[i].mass));
          }
          bodies[i].resetForce();
        }

        // if theres collisions, then make new body and remove the collided bodies
        checkCollision(bodies, collisions, collisionPairCount);
        // increment the current time with the time step.
        // write to csv the number of bodies at this state if it's changed.
        if ((BeforeBodyCount - NumberOfBodies) != 0){
          countWrite(std::to_string(t) + "," + std::to_string(NumberOfBodies) + "\n");
        }
      }
    } 
   t += currentTimestep;
  }

  // check if the arguments are valid
  int verifyArguments(int argc, char** argv){
    if (argc==1) {
      std::cout << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
                << std::endl
                << "Examples:" << std::endl
                << "100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
                << "100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
                << "100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
                << std::endl
                << "In this naive code, only the first body moves" << std::endl;

      return -1;
    }
    else if ( (argc-2)%7!=0 ) {
      std::cout << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
      return -2;
    }
    return 0;
  }

  // Saves running time
  void timeHandler(){
    clock_t tEnd = clock();

    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    std::cout << now->tm_mday << '/'
        << (now->tm_mon + 1) << '/'
        << (now->tm_year + 1900) << " "
        << std::setw(2) << std::setfill('0')
        << now->tm_hour << ":"
        << std::setw(2) << std::setfill('0')
        << now->tm_min << ":"
        << std::setw(2) << std::setfill('0')
        << now->tm_sec
        << std::endl;

    double timeTaken = (tEnd - tStart) * 0.000001;

    // load file to write to
    std::ofstream outfile;

    std::cout << numberOfCores << " Time taken: " << timeTaken << "s\n";

    // save to file
    outfile.open("timerResults.txt", std::ios_base::app);
    outfile << now->tm_mday << '/'
        << (now->tm_mon + 1) << '/'
        << (now->tm_year + 1900) << " "
        << std::setw(2) << std::setfill('0')
        << now->tm_hour << ":"
        << std::setw(2) << std::setfill('0')
        << now->tm_min << ":"
        << std::setw(2) << std::setfill('0')
        << now->tm_sec
        << " - "
        << "#Bodies: " << initialNumberOfBodies
        << ", Timestep: " << defaultTimeStepSize
        << ", Runtime: " << timeTaken << "s\n";
  }

  // Starts the space body simulations.
  int performSpaceBodies(Body* bodies){
    initialNumberOfBodies = NumberOfBodies;
    timeStepsSinceLastPlot = 0;
    const int plotEveryKthStep = 100;
   
    while (t<=tFinal) {
      updateBodies(bodies);
      timeStepsSinceLastPlot++;
      // this is used to watch if the debug is called to check for collisions.
      if (isCollided && isCsvCollisionWrite){
        break;
      }
      if ((timeStepsSinceLastPlot%plotEveryKthStep==0) && useParaview) {
        printParaviewSnapshot(timeStepsSinceLastPlot/plotEveryKthStep, bodies);
      }
    }
    timeHandler();
    return 0;
  }

// INITIALISE FUNCTIONS -------------------------

  // Generates a random set of bodies.
  void generateRandomBodies(Body* b){

    // I dont understand why so many people upvoted this answer. It is mathematically incorrect. RAND_MAX is a very small number (typically 2^16). That means that from 23 bits of the floating point you make only 15 random. The others will be probably zero. You will indeed get random numbers in uniform distribution but of low precision. For example your random generator can generate 0.00001 and 0.00002 but cannot generate 0.000017. So you have a uniform distribution but of low precision (256 times less precision than the actual floating point). – DanielHsH Aug 14 '14 at 8:25

    double maxFloat = 1.0;

    for (int i = 0; i < NumberOfBodies; i++){
      double random_value[7];
      for(int j = 0; j < 7; j++){
        // generate random double
          double f = (double)rand() / RAND_MAX;
          f = fMin + f * (fMax - fMin);

        // spread range from -1 to 1 non inclusive and add to array
        random_value[j] = f;
      }
      // initialise the bodies values.
      for(int k = 0; k < 3; k++){
        b[i].x[k] = round( random_value[k] * 1000.0 ) / 1000.0; // round position to 3dp
        b[i].v[k] = round( random_value[k+3] * 1000.0 ) / 1000.0; // round velocity to 3dp
        b[i].force[k] = 0;
      }
      // make sure the mass has a positive value!
      b[i].mass = fabs(random_value[6]);
    }
  }

  // Runs the random bodies simulation
  void runRandomBodies(){
    // test with random bodies
      NumberOfBodies = NumberOfRandomBodies;
      tFinal = randomSimTFinal;
      Body bodies[NumberOfBodies];

      if (useParaview){
        openParaviewVideoFile();
        printParaviewSnapshot(0, bodies);
      }
      
      generateRandomBodies(bodies);
      performSpaceBodies(bodies);

      if (useParaview){closeParaviewVideoFile();}
  }

  // Initiate runtime
  int main(int argc, char** argv) {
    // set random seed (should be different every time!)
    srand(seed);

    // check for argument size
    if (argc > 1){
      RandomBodies = false;
    }
    if (RandomBodies == false){
      if (verifyArguments(argc,argv) == 0){
      // initiate bodies array
      NumberOfBodies = (argc-2) / 7; // count the number of bodies.
      Body bodies[NumberOfBodies];
      // initial setup.
      
      int loopCap = 1;
      if (isCsvCollisionWrite && collisionIterate){
        loopCap = numberOfIterations;
      }

      if (isCsvCollisionWrite){
        // write header for csv if enabled
        collisionWrite("Timestep,");
        collisionWrite("distance\n");
      }

      std::cout << "Performing Space Bodies" << std::endl;
      countWrite("Time, Number of Bodies\n");
      countWrite(std::to_string(t) + "," + std::to_string(NumberOfBodies) + "\n");
      
      for (int s = 0; s < loopCap; s++){
        // reset values
        for (int i = 0; i < NumberOfBodies; i++){bodies[i].reset();}
        NumberOfBodies = (argc-2) / 7; // count the number of bodies.
        t = 0;
        tFinal = 0;
        isCollided = false;

        // set up values again
        setUp(argc,argv,bodies);
        
        if (useParaview){
          openParaviewVideoFile();
          printParaviewSnapshot(0, bodies);
        }

        // perform space loops.
        performSpaceBodies(bodies);
        defaultTimeStepSize = defaultTimeStepSize / 2;
      }

      if (useParaview){closeParaviewVideoFile();}
      }
    } else {
      runRandomBodies();
    }
    // Close the files that were opened in memory.
    if (isCsvBodyCountWrite){csvBodyCountFile.close();}
    if (isCsvCollisionWrite){csvCollisionFile.close();}
    return 0;
  }
