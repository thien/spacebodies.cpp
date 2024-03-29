// Translate this file with
//
// g++ -O3 --std=c++11 spacebodies.cpp -o spacebodies
// or in parallel
// g++ -O3 --std=c++11 -fopenmp spacebodies.cpp -o spacebodies
//
// Run it with
//
// ./spacebodies
// 
// if you're planning to utilise the random bodies options (configure below)
// 
// ./spacebodies -r
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

// CONFIGURATION OPTIONS

bool adaptiveTimeStepCheck = true; // If true, then use adaptive timestep

// If you're using paraview,
// make sure that there is a 'paraview' folder in the same directory as this file!
bool useParaview = true; // if true, write paraview files.
bool saveEveryStepToParaview = false; // if true, saves every step. otherwise it saves contingent on equal increments of time.

double defaultTimeStepSize = 0.00001;
double smallSizeLimit = 1e-8; // If variables are smaller than this then it might as well be zero.

int numberOfIterations = 20; // When using the collision iteration (for different timesteps)
bool isCsvCollisionWrite = false; // if set to true, it will generate a csv of two bodies and data to show their collision.
bool terminateOnCollision = false; // if set to true, terminate the program after a collision
bool collisionIterate = false; // iterates through multiple rounds, halving the timestep size as it goes.
bool writeTimerCount = false; //if true, saves runtime logs of every simulation

bool isCsvBodyCountWrite = false; // if true, will write a csv that counts the number of bodies over time.

// writer variables
bool printBodiesInfo = false; // prints the bodies and their data
bool printTimestampInfo = true;// print timestamp if true

double referenceDistance = -1.95; // If measuring error, use this value as the reference!

// Tools to manage random spacebodies
// random bodies are called with the -r flag..
// i.e ./spacebodies -r

int NumberOfRandomBodies = 100;
double finalRandomSimulationTime = 10.0;
double fMin = 1; // random double min
double fMax = 2; // random double max
bool bodiesHaveMass = true; // if false, then the bodies have negligible mass

// Seed value (used to generate random values for the bodies.)
double seed = 12321321;








// SETUP VARIABLES-------------------------------
  // DON'T EDIT THESE..

  bool RandomBodies = false;
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
  time_t startTime = time(0);   // get time now

  // timer variables
  clock_t genRandomBodiesStart = clock();
  clock_t genRandomBodiesEnd = clock();
  clock_t initEnd = clock();
  clock_t parallelEndTime = clock();
  clock_t taskEndTime = clock();  


  std::string startTimeString(){
    struct tm * now = localtime( & startTime );
    
    std::stringstream date_string;
    date_string << now->tm_mday << '-'
        << (now->tm_mon + 1) << '-'
        << (now->tm_year + 1900) << "_"
        << std::setw(2) << std::setfill('0')
        << now->tm_hour << ":"
        << std::setw(2) << std::setfill('0')
        << now->tm_min << ":"
        << std::setw(2) << std::setfill('0')
        << now->tm_sec;
    // std::cout << date_string.str();
    return date_string.str();
  }

// CSV WRITER HANDLERS --------------------------

  std::ofstream csvBodyCountFile;
  std::ofstream csvCollisionFile;

  // helper functions to write to csv.
  void countWrite(std::string text){
    if (isCsvBodyCountWrite){
      csvBodyCountFile.open("bodycount.csv", std::ios_base::app);
      csvBodyCountFile << text;
      csvBodyCountFile.close();
    }
  }
  void collisionWrite(std::string text){
    if (isCsvCollisionWrite){
      csvCollisionFile.open("collision.csv", std::ios_base::app);
      csvCollisionFile << text;
      csvCollisionFile.close();
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

      void resetForce(){
        for (int i = 0; i < 3; i++) {force[i] = 0;}
      }

      void print(int bodyID){
        // prints body statistics.
        printf("Body %4d: %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f \n", bodyID, x[0], x[1], x[2], v[0], v[1], v[2], mass);
      }
  };

// PARAVIEW FUNCTIONS ---------------------------

  /*
  * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
  */

  void printParaviewSnapshot(int counter, Body * bodies) {

    std::stringstream filename;
    filename << startTimeString() << "_result-" << counter <<  ".vtp";
    std::stringstream filepath;
    filepath << "paraview/" << filename.str();
    std::ofstream out( filepath.str().c_str() );
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << "\t<Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "\t\t<Points>" << std::endl
        << "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
        << std::endl;

    for (int i=0; i<NumberOfBodies; i++) {
      out << "\t\t\t\t"
          << bodies[i].x[0]
          << " "
          << bodies[i].x[1]
          << " "
          << bodies[i].x[2]
          << std::endl;
    }

    out << "\t\t\t</DataArray>" << std::endl
        << "\t\t</Points>" << std::endl
        << "\t</Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>"  << std::endl;

    videoFile << "\t\t<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
  }

  void openParaviewVideoFile() {
    std::stringstream pvdFile;
    pvdFile << "paraview/" << startTimeString() << "_result.pvd";

    videoFile.open( pvdFile.str() );
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "\t<Collection>"
              << std::endl;
  }

  void closeParaviewVideoFile() {
    videoFile << "\t</Collection>" << std::endl
              << "</VTKFile>" << std::endl;
  }

// COLLISION HANDLERS ---------------------------

  void collisionDebug(Body a, Body b){
    if (isCsvCollisionWrite){
      double pos = a.x[0]-referenceDistance;
      
      std::ostringstream strs;
      strs << pos;
      std::string posString = strs.str();

      std::ostringstream str2;
      str2 << currentTimestep;
      std::string ts = str2.str();

      std::ostringstream str3;
      str3 << currentTimestep / pos;
      std::string curTS_pos_Str = str3.str();

      std::ostringstream str4;
      str4 << NumberOfBodies;
      std::string bodyCountStr = str4.str();

      if (isCsvCollisionWrite){
        std::cout << "COLLISION @ t="<< t<<", Timestep: "<< currentTimestep << ": Error: " << pos << ", Location: " << a.x[0] << "; " << b.x[0] << ", Ratio: " << (currentTimestep/pos) <<", Steps: "<< timeStepsSinceLastPlot << ", Range: " << abs(b.x[0] - a.x[0]) << std::endl;
      }

      std::cout <<currentTimestep<<"," << pos << "," << a.x[0] << "," << b.x[0] << "," << (currentTimestep/pos) <<","<< timeStepsSinceLastPlot << "," << abs(b.x[0] - a.x[0]) << std::endl;
      
      collisionWrite(ts+ ", " + posString + ", "+curTS_pos_Str+"," + bodyCountStr + "\n");
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
    }
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
            // we want the new distance to be around 1/5 of the distance they will cover.
            badDistanceRatio = (((newDist-currentDistance) / currentDistance) < 0.2);
          }
        }
      }
    }
    return timestep;
  }

  // function that updates the positions of the particles in space
  void updateBodies(Body* bodies) {
    int i,j,k; //incrementer variables
    currentTimestep = defaultTimeStepSize;
    if (printBodiesInfo){printf("\n\n");}
    if (printTimestampInfo){
      if (printBodiesInfo){
        // if we're printing bodies too then no point using carriage returns
        // since C++ doesn't support multiple line manipulation
        printf ("T: %012.10f, # Bodies: %1.0d \n", t, NumberOfBodies);
      } else {
        // inline print to save lines on the terminal
         std::cout << "# Bodies: " << NumberOfBodies << ", Time: " << std::fixed << std::setprecision(9) << t << "\r"<< std::flush;
      }
    }


    // initiate positions of the shortest body positions
    double closestDistance = 999999;
    int closestPair1 = 0;
    int closestPair2 = 0;
    // initiate values we'll use to determine collided bodies
    int collisions[NumberOfBodies][2];
    int collisionPairCount = 0;
    int BeforeBodyCount = NumberOfBodies;

    #pragma omp parallel default(none) private(j,k) shared(i, parallelEndTime, isCollided, smallSizeLimit, adaptiveTimeStepCheck, BeforeBodyCount, printBodiesInfo, t, bodies, currentTimestep, closestDistance, closestPair1, closestPair2, collisions, collisionPairCount, NumberOfBodies)
    {
      
      // Calculate distances and forces for each body
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
      parallelEndTime = clock();
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
   taskEndTime = clock();
  }

  // Saves running time to file
  void timeHandler(){
    clock_t tEnd = clock();

    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    // std::cout << now->tm_mday << '/'
    //     << (now->tm_mon + 1) << '/'
    //     << (now->tm_year + 1900) << " "
    //     << std::setw(2) << std::setfill('0')
    //     << now->tm_hour << ":"
    //     << std::setw(2) << std::setfill('0')
    //     << now->tm_min << ":"
    //     << std::setw(2) << std::setfill('0')
    //     << now->tm_sec
    //     << std::endl;


    if (writeTimerCount){
      double timeMult = 0.000001;
      double timeTaken = (tEnd - tStart) * timeMult;

      double randomBodyGenTime = (genRandomBodiesEnd - genRandomBodiesStart) * timeMult;
      double initSetupTime = (initEnd - tStart) * timeMult;
      double parallelTime = (parallelEndTime - initEnd) * timeMult;
      double SerialPostTime = (taskEndTime - parallelEndTime) * timeMult;

      std::cout <<"Init Setup Time: " << initSetupTime;
      std::cout << "s, BodyGen Time: " << randomBodyGenTime;
      std::cout << "s, Parallel Time: " << parallelTime;
      std::cout << "s, SerialPost Time: " << SerialPostTime;
      std::cout << "s\n";
      std::cout << "Time taken: " << timeTaken << "s\n";

      // load file to write to
      std::ofstream outfile;

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
          << ", CPU Time: " << timeTaken << "s\n";
    }
  }

  // Starts the space body simulations.
  int performSpaceBodies(Body* bodies){
    // add start time.
    initEnd = clock();
    initialNumberOfBodies = NumberOfBodies;
    timeStepsSinceLastPlot = 0;
    const int plotEveryKthStep = 100;
   
    // repeat the loop for every step until we reach the end time.
    while (t<=tFinal) {
      updateBodies(bodies);
      timeStepsSinceLastPlot++;
      // this is used to watch if the debug is called to check for collisions.
      if (isCollided && isCsvCollisionWrite && terminateOnCollision){
        break;
      }
      // update paraview files if enabled.
      if (useParaview){
        if(saveEveryStepToParaview){
          printParaviewSnapshot(timeStepsSinceLastPlot, bodies);
        } else {
          if (timeStepsSinceLastPlot%plotEveryKthStep==0){
            printParaviewSnapshot(timeStepsSinceLastPlot/plotEveryKthStep, bodies);
          }
        }
      }
    }
    // end time is handled here.
    timeHandler();
    return 0;
  }

// INITIALISE FUNCTIONS -------------------------

  // Generates a random set of bodies.
  void generateRandomBodies(Body* b){
    for (int i = 0; i < NumberOfBodies; i++){
      double random_value[7];
      for(int j = 0; j < 7; j++){
        // generate random double
        double f = (double)rand() / RAND_MAX;
        f = fMin + f * (fMax - fMin);
        random_value[j] = f;
      }
      // initialise the bodies values.
      for(int k = 0; k < 3; k++){
        b[i].x[k] = random_value[k];
        b[i].v[k] = random_value[k+3] * 100;
        b[i].force[k] = 0;
      }
      if (bodiesHaveMass){
        b[i].mass = fabs(random_value[6]); // make sure the mass has a positive value!
      } else {
        b[i].mass = 0.0000000001;
      }
    }
  }

  // check validity of arguments
  int verifyArguments(int argc, char** argv){

    if (argc==1) {
      std::cout << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
                << std::endl
                << "Examples:" << std::endl
                << "100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
                << "100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
                << "100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
                << std::endl;
      return -1;
    } else {
      for (int i = 0; i < argc; i++){
        // look at flags.
        std::string arg = argv[i];
        if (arg == "-r"){
          std::cout << "Random Body Flag is found!" << std::endl;
          RandomBodies = true;
        }
      }

      if ( ((argc-2)%7!=0) && (RandomBodies == false)) {
        std::cout << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
        return -2;
      }
    }
    return 0;
  }

  // Initiate runtime
  int main(int argc, char** argv) {
    srand(seed);
    int loopCap = 1;
    int defaultBodyCount;

    // make sure the arguments are valid prior to running.
    // Also look at the arguments prior.
    if (verifyArguments(argc,argv) == 0){
      std::cout << "Getting Body Count.. ";

      if (RandomBodies){
        defaultBodyCount = NumberOfRandomBodies;
      } else {
        defaultBodyCount = (argc-2)/7;
      }

      std::cout << "Done; Found " << defaultBodyCount << " bodies." << std::endl;
      // assign default body count to NumberOfBodies variable (used everywhere except main)
      NumberOfBodies = defaultBodyCount;
      // initiate body array
      Body bodies[NumberOfBodies];
    
      if (isCsvCollisionWrite){
        // write header for csv if enabled
        collisionWrite("Timestep, distance\n");
        if (collisionIterate){loopCap = numberOfIterations;}
      }

      // print to countWriter;
      countWrite("Time, Number of Bodies\n");
      countWrite(std::to_string(t) + "," + std::to_string(defaultBodyCount) + "\n");
      
      // initate run.
      for (int s = 0; s < loopCap; s++){
        std::cout << "Resetting values for next iteration if needed.." << std::endl;
        for (int i = 0; i < defaultBodyCount; i++){bodies[i].reset();}
      
        // generate bodies.
        if (RandomBodies){
          std::cout << "Initiating Random Bodies.. ";
          genRandomBodiesStart = clock();
          generateRandomBodies(bodies);
          genRandomBodiesEnd = clock();
          tFinal = finalRandomSimulationTime;
        } else {
          std::cout << "Initiating argument bodies.. ";
          setUp(argc,argv,bodies);
        }
        std::cout << " Done." << std::endl;
        
        if (useParaview){
          std::cout << "Initiating Paraview files.." << std::endl;
          openParaviewVideoFile();
          printParaviewSnapshot(0, bodies);
        }

        std::cout << "Starting The Space Bodies Program.." << std::endl;
        performSpaceBodies(bodies);

        if (loopCap > 1){
          // half the timestep for the next iteration (if there is one.)
          defaultTimeStepSize = defaultTimeStepSize / 2;
          t = 0;
          isCollided = false;
        }
      }

      if (useParaview){
        std::cout << "Closing paraview files.." << std::endl;
        closeParaviewVideoFile();
      }
    }
  
    // Close the files that were opened in memory.
    return 0;
  }
