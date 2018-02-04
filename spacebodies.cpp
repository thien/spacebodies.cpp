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
#include <iostream>
#include <string>
#include <cstdlib>
#include <math.h>
#include <unordered_map>
// library is used to deduce leading zeros print
#include <iomanip>


double t = 0;
double tFinal = 0;
bool adaptiveTimeStepCheck = true;
double defaultTimeStepSize = 0.01;
int NumberOfBodies = 0;
double smallSizeLimit = 1e-8;

// ---------------------------

// Body Class; It'll be useful for when I need to make adjustments.

class Body {
  public:
    double force[3]; //body force
    double x[3];   // displacement of body
    double v[3];   // velocity of body
    double mass;    // body mass
    double nx[3];  // future displacement of body
    double nv[3];  // future velocity of body
    double nm;   // future body mass

    void assignNewValues(){
      // assigns future values to be current values.
      for (int i = 0; i < 3; i++) {
          x[i] = nx[i];
          v[i] = nv[i];
      }
      mass = nm;
      // also reset the force values when we're done.
      resetForce();
    }

    // If the bodies force is small enough then we set it to zero.
    void capForceLimit(){
      for (int k=0; k<3; k++){
        if (force[k] < smallSizeLimit){
          force[k] = 0;
        }
      }
    }

    void printForce(){
      std::cerr << "    Force:  "<<force[0]<< ", "<<force[1]<< ", "<<force[2]<<std::endl;
    }

    void resetForce(){
      // printf("resetting forces \n");
      for (int i = 0; i < 3; i++) {
          force[i] = 0;
      }
    }

    void print(int bodyID){
      // literally prints body statistics.
      // Values are (ID, x,y,z, v(x), v(y), v(z), mass)
      printf("Body %4d: %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f  %+010.6f \n",
        bodyID, x[0], x[1], x[2], v[0], v[1], v[2], mass);
    }
};

// ---------------------------

// Floating Point Class;
  // Useful for self learning (+ prep for exams
  const int maxS = 10000;

  class myFloat {
    public:
      int s,e;
      myFloat normalise(myFloat f){
        myFloat result;
        result.s = f.s;
        result.e = f.e;
        while(result.s > maxS){
          result.s = result.s / 10;
          result.e = result.e + 1;
        }
        while(result.s < maxS/10){
          result.s = result.s * 10;
          result.e = result.e - 1;
        }
      return result;
      }
      myFloat add(myFloat a, myFloat b){
        while (a.e < b.e){
          b.e = b.e-1;
          b.s = b.s*10;
        }
        myFloat result;
        result.s = a.s + b.s;
        result.e = a.e;
        return normalise(result);
      }
  };


// --------------------------- TAYLOR STUFF

  // // Taylor Series (currently Sine as example)
  // // convert to radians to degree
  // double toRadians(double angdeg){
  //   //x is in radians
  //   const double PI = 3.14159265358979323846;
  //   return angdeg / 180.0 * PI;
  // }
  // //factorial function
  // double fact(double x){
  //   // calculate factorial for denominator
  //   if (x==0 || x==1) {
  //     return 1;
  //   } else {
  //     return (x * fact(x - 1));
  //   }
  // }
  // //mySin function
  // double mySin(double x){
  //   double sum = 0.0;
  //   for(int i = 0; i < 9; i++){
  //     double top = pow(-1, i) * pow(x, 2 * i + 1);  //calculation for nominator
  //     double bottom = fact(2 * i + 1);              //calculation for denominator
  //     sum = sum + top / bottom;                     //1 - x^2/2! + x^4/4! - x^6/6!
  //   }
  //   return sum;
  // }
  // // Test Code for sine
  // void runSin(){
  //   double param = 45, result;
  //   result= mySin(toRadians(param));
  // }

// ---------------------------
/*
 * Paraview Functions
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */

  std::ofstream videoFile;

  void printParaviewSnapshot(int counter, Body * bodies) {
    std::stringstream filename;
    filename << "result-" << counter <<  ".vtp";
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
          << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>"  << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
  }

  void openParaviewVideoFile() {
    videoFile.open( "result.pvd" );
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "<Collection>";
  }

  void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
  }

// -------------------------------------------

// Function to join two collided bodies and joins them together to create
// A merged body.
Body joinBodies(Body x, Body y){
  // joins the forces and the masses of two bodies.
  // The new body shall have a new mass that is equal to the sum of the masses of the original bodies and an averaged velocity (derive the formulae and present it in the report; you might want to implement some momentum-preserving scheme, i.e. mass-averaged result velocities, but the assignment is about the coding not about realistic physics, so a simple average of velocity components does the job).

  Body z;
  // combine mass.
  z.mass = x.mass + y.mass;
  for (int i; i < 3; i++){
    // get position of one of the colluded bodies
    z.x[i] = x.x[i];
    // calculate velocity of new item. (Mass averaged velocities)
    // std::cerr << "V"<<i<<": "<< z.v[i] << std::endl;
    z.v[i] = (x.v[i]/x.mass) + (y.v[i]/y.mass);
    // std::cerr << "V"<<i<<": "<< z.v[i] << std::endl;
  }
  return z;
}

void checkCollision(Body b[]){
  // create space so we can find out which bodies are colliding with each other
  int positions[NumberOfBodies][2];
  int collisionPairCount = 0;
  std::unordered_map<std::string,bool> used_bodies;

  // find all possible collisions between particles.
  for (int i = 0; i < NumberOfBodies; i++){
    for (int j = 0; j < i; j++){
      if (i != j){ // don't calculate distance from itself to itself.
        double pos = sqrt(
          (b[j].x[0]-b[i].x[0]) * (b[j].x[0]-b[i].x[0]) +
          (b[j].x[1]-b[i].x[1]) * (b[j].x[1]-b[i].x[1]) +
          (b[j].x[2]-b[i].x[2]) * (b[j].x[2]-b[i].x[2])
        );
        if (pos <= smallSizeLimit){
          // Collision means the bodies are closer than 1e-8.
          positions[collisionPairCount][0] = i;
          positions[collisionPairCount][1] = j;
          collisionPairCount++;
        }
      }
    }
  }

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

// -------------------------------------------

// Core Functions

// DONE
int calculateNumberOfBodies(int argc){
  // this represents the initial number of bodies in space.
  return (argc-2) / 7;
}

// DONE
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
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
    printf("Body %5.0d \t %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f \n", i, b[i].x[0], b[i].x[1], b[i].x[2], b[i].v[0], b[i].v[1], b[i].v[2], b[i].mass);
  }
  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
}

// TODO
double updateTimeStep(double beforeTS, Body a, Body b){

  double currentDistance = sqrt(
    (b.x[0]-a.x[0]) * (b.x[0]-a.x[0]) +
    (b.x[1]-a.x[1]) * (b.x[1]-a.x[1]) +
    (b.x[2]-a.x[2]) * (b.x[2]-a.x[2])
  );

  double timestep = beforeTS;
  
  // check if that timestep is tiny.

  // we have a distance that is quite small.
  // we need to see whether the new distance will collide with them or go past them
  // get the timestep and guess a new distance that would be at least 1/4 and at most 1/2 the current distance between them.

  // is the timestep too small? (your cap is 1e^-10)


  // if (currentDistance < 1){
  //   double newAX[3];
  //   double newBX[3];
  //   double newDist;
  //   bool notWithinRegion = true;


  //   for (int i = 0; i < 3; i++){
  //     newAX[i] = (0.0000001 * a.v[i]) + a.x[i];
  //     newBX[i] = (0.0000001 * b.v[i]) + b.x[i];
  //   }
  //   newDist = sqrt(
  //     (newBX[0]-newAX[0]) * (newBX[0]-newAX[0]) +
  //     (newBX[1]-newAX[1]) * (newBX[1]-newAX[1]) +
  //     (newBX[2]-newAX[2]) * (newBX[2]-newAX[2])
  //   );
  //   if(newDist - currentDistance > 0){
  //     // if they're going further away then we don't need to bother about it.
  //     notWithinRegion = false;
  //   }

  //   // if they are within the region:
  //   while (notWithinRegion){
  //     timestep = timestep / 2;
  //     std::cerr << " Timestep halved:  "<<timestep<< std::endl;
  //     for (int i = 0; i < 3; i++){
  //       newAX[i] = (timestep * a.v[i]) + a.x[i];
  //       newBX[i] = (timestep * b.v[i]) + b.x[i];
  //     }
  //     newDist = sqrt(
  //       (newBX[0]-newAX[0]) * (newBX[0]-newAX[0]) +
  //       (newBX[1]-newAX[1]) * (newBX[1]-newAX[1]) +
  //       (newBX[2]-newAX[2]) * (newBX[2]-newAX[2])
  //     );
  //     std::cerr << " New Dist:  "<<newDist << ", " << currentDistance << std::endl;
      

  //     notWithinRegion = ((newDist / currentDistance) > 0.3);
  //     notWithinRegion = ((newDist / currentDistance) < 0.4);

  //     // notWithinRegion = (newDist * 2 <= currentDistance);
  //   }

    // get the current value and see if we can get it to half the distance
  // }
  if (timestep >= smallSizeLimit){

      double EPS = 1e-8;
      for (int i = 0; i < 3; i++){

        if (a.v[i] > 0 and (abs(timestep * a.force[0]/a.mass)/a.v[0] > EPS)){
          timestep = timestep/2;
        } 
      }
  }
  return timestep;
}

// function that updates the positions of the particles in space
void updateBodies(Body* bodies) {
  printf ("\n\nTime: %4.8f, NumberOfBodies: %1.0d \n", t, NumberOfBodies);
  double timestep = defaultTimeStepSize;

  // initiate positions of the shortest body positions
  double closestDistance = 999999;
  int closestIndex1 = 0;
  int closestIndex2 = 0;

  for (int j=0; j<NumberOfBodies; j++) {
    // now we can print the position of the space body
    bodies[j].print(j);
    // only calculate half the calculations!
    for (int i=0; i<j; i++) {
      if (i != j){ // make sure it doesn't interact with itself.
        // calculate distance
        double distance = sqrt(
          (bodies[j].x[0]-bodies[i].x[0]) * (bodies[j].x[0]-bodies[i].x[0]) +
          (bodies[j].x[1]-bodies[i].x[1]) * (bodies[j].x[1]-bodies[i].x[1]) +
          (bodies[j].x[2]-bodies[i].x[2]) * (bodies[j].x[2]-bodies[i].x[2])
        );
        // check if this is the closest body to date
        if (closestDistance > distance){
          closestDistance = distance;
          closestIndex1 = i;
          closestIndex2 = j;
        }
        // calculate combined mass
        double combinedMass = bodies[i].mass * bodies[j].mass;
        // update force values (for x,y,z)
        double calc5 = distance/distance/distance;
        double calc4 = calc5/combinedMass;
        // std::cerr <<"D: "<<distance << " C1:  "<<calc1 <<", C2: "<<calc2 <<", C3: "<<calc3 <<", C4: "<<calc4 <<", C5: "<<calc5 << std::endl;
        for (int k = 0; k < 3; k++){
          bodies[j].force[k] += abs(bodies[i].x[k]-bodies[j].x[k]) * calc4;
          bodies[i].force[k] += abs(bodies[i].x[k]-bodies[j].x[k]) * calc4;
        }
      }
    }
    // check if the force is too small; if it is set it to 0 to make life easier.
    bodies[j].capForceLimit();
  }

  if (adaptiveTimeStepCheck == true){
    // update timestep if needed
    timestep = updateTimeStep(timestep, bodies[closestIndex2], bodies[closestIndex1]);
  }

  // once we've also calculated the timestep
  for (int j=0; j<NumberOfBodies; j++){
    // update position and velocity of body
    for (int k=0; k<3; k++){
      // update the position and velocity of the body.
      bodies[j].nx[k] = bodies[j].x[k] + (timestep * bodies[j].v[k]);
      bodies[j].nv[k] = bodies[j].v[k] + (timestep * (bodies[j].force[k] / bodies[j].mass));
    }
    bodies[j].nm = bodies[j].mass;
  }

  for (int j=0; j<NumberOfBodies; j++){
    // make every space body assign their new values.
    bodies[j].assignNewValues();
  }
  // check for any collisions
  // if theres collisions, then make new body and remove the collided bodies
  checkCollision(bodies);
  // increment the current time with the time step.
  t += timestep;
  // std::cout << "\x1B[2J\x1B[H";
}

// DONE
int verifyArguments(int argc, char** argv){
  if (argc==1) {
    std::cerr << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
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
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }
  return 0;
}

// Starts the space body simulatiions.
int performSpaceBodies(Body* bodies){
  std::cerr << "Performing Space Bodies" << std::endl;

  int timeStepsSinceLastPlot = 0;
  const int plotEveryKthStep = 100;
  while (t<=tFinal) {
    updateBodies(bodies);
    timeStepsSinceLastPlot++;
    // if (timeStepsSinceLastPlot%plotEveryKthStep==0) {
    //   printParaviewSnapshot(timeStepsSinceLastPlot/plotEveryKthStep);
    // }
  }
  return 0;
}

// -------------------------------------------
// Generation of Space Bodies
// -------------------------------------------

// Generates a random set of bodies.
void generateRandomBodies(Body* b){

  // I dont understand why so many people upvoted this answer. It is mathematically incorrect. RAND_MAX is a very small number (typically 2^16). That means that from 23 bits of the floating point you make only 15 random. The others will be probably zero. You will indeed get random numbers in uniform distribution but of low precision. For example your random generator can generate 0.00001 and 0.00002 but cannot generate 0.000017. So you have a uniform distribution but of low precision (256 times less precision than the actual floating point). – DanielHsH Aug 14 '14 at 8:25

  float maxFloat = 1.0;

  for (int i = 0; i < NumberOfBodies; i++){
    float random_value[7];
    for(int j = 0; j < 7; j++){
      // generate random value between 0-1 (as a float)
      float random_val = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      // spread range from -1 to 1 non inclusive and add to array
      random_value[j] = (random_val - 0.5) * 2;
    }
    // assign the random values to the body
    for(int k = 0; k < 3; k++){
      b[i].x[k] = random_value[k];
      b[i].v[k] = random_value[k+3];
    }
    // make sure the mass has a positive value!
    b[i].mass = fabs(random_value[6]);
  }
}

// Generates a set of bodies that collide with each other.
void generateCollidingBodies(Body* b){
  b[0].x[0] = 1;
  b[0].x[1] = 1;
  b[0].x[2] = 1;
  b[0].v[0] = 0;
  b[0].v[1] = 0;
  b[0].v[2] = 0;
  b[0].mass = 1;

  b[1].x[0] = 1;
  b[1].x[1] = 1;
  b[1].x[2] = 1;
  b[1].v[0] = 0;
  b[1].v[1] = 0.00001;
  b[1].v[2] = 0;
  b[1].mass = 1;
}

// -------------------------------------------

// Runs the colliding bodies simulation
void runCollidingBodies(){
  // test with random bodies
    NumberOfBodies = 2;
    tFinal = 0.01;
    Body bodies[NumberOfBodies];
    generateCollidingBodies(bodies);
    performSpaceBodies(bodies);
}

// Runs the random bodies simulation
void runRandomBodies(){
  // test with random bodies
    NumberOfBodies = 40;
    tFinal = 20.0;
    Body bodies[NumberOfBodies];
    generateRandomBodies(bodies);
    performSpaceBodies(bodies);
}

// Initiate runtime
int main(int argc, char** argv) {
  // check for arguments size.
  bool RandomBodies = true;
  if (argc > 1){
    RandomBodies = false;
  }
  if (RandomBodies == false){
    if (verifyArguments(argc,argv) == 0){
    // initiate bodies array
    NumberOfBodies = calculateNumberOfBodies(argc);
    Body bodies[NumberOfBodies];
    // initial setup.
    setUp(argc,argv,bodies);

    // openParaviewVideoFile();
    // printParaviewSnapshot(0);
    // perform space loops.
    performSpaceBodies(bodies);
    // close the paraview files.
    // closeParaviewVideoFile();
    // exit safely.
    }
  } else {
    runCollidingBodies();
  }
  return 0;
}
