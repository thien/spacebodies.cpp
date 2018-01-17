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

double t = 0;
double tFinal = 0;
bool adaptiveTimeStepCheck = true;
double defaultTimeStepSize = 1;
double smallestDiffCoefficent = 100;
int NumberOfBodies = 0;

// ---------------------------

// Body Class; It'll be useful for when I need to make adjustments.

class Body {
  
  public:
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


// ---------------------------

// Helper Functions (Because they're not included in C++)

// Binary search for integer in list of integers.
int BinarySearch(int* arr, int n, int search){
  int first, last, middle, location;
  // default location value to -1   
  location = -1;
  // init holder values
  first = 0;
  last = n-1;
  middle = (first+last)/2;

  while (first <= last){
    if(arr[middle] < search){
      first = middle + 1;
    }
    else if(arr[middle] == search){
      // we found the location
      location = middle + 1;
      break;
    }
    else {
      last = middle - 1;
    }
    middle = (first + last)/2;
  }
  if(first > last){
    // not found.
  }
  return location;

  // Example runtime:
    // int arr[50];//init array size
    // int n = 20; // refers to the number of items in the array
    // for (int i=0; i<n; i++){
    //     arr[i] = i*2;
    // }
    // int search = 4; //search for this value
    // int k = BinarySearch(arr, n, search);
    // printf("%10d", k);
}

// Boolean model of Binary Search
bool ItemInArray(int* arr, int n, int search){
  if (BinarySearch(arr,n,search) != -1){
    return false;
  } else {
    return true;
  }
}

// TODO
void removeObjectsFromArray(){
  // makes new arrays for mass, x, and v
  // copy pointers from previous mass,x,v
  // remove pointers that need to be removed
  // add new pointer of object.
}

// ---------------------------

// Taylor Series (currently Sine as example)
// convert to radians to degree
double toRadians(double angdeg){                                       
  //x is in radians
  const double PI = 3.14159265358979323846;
  return angdeg / 180.0 * PI;
}
//factorial function 
double fact(double x){
  // calculate factorial for denominator
  if (x==0 || x==1) {
    return 1;
  } else {
    return (x * fact(x - 1));
  }
}
//mySin function
double mySin(double x){
  double sum = 0.0;
  for(int i = 0; i < 9; i++){
    double top = pow(-1, i) * pow(x, 2 * i + 1);  //calculation for nominator
    double bottom = fact(2 * i + 1);              //calculation for denominator
    sum = sum + top / bottom;                     //1 - x^2/2! + x^4/4! - x^6/6!
  }
  return sum;
}
// Test Code for sine
void runSin(){
  double param = 45, result;
  result= mySin(toRadians(param));
}

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

// ---------------------------

// DONE
Body joinBodies(Body x, Body y){
  // joins the forces and the masses of two bodies.
  // The new body shall have a new mass that is equal to the sum of the masses of the original bodies and an averaged velocity (derive the formulae and present it in the report; you might want to implement some momentum-preserving scheme, i.e. mass-averaged result velocities, but the assignment is about the coding not about realistic physics, so a simple average of velocity components does the job).

  Body z;
  z.mass = x.mass + y.mass;
  for (int i; i < 3; i++){
    z.x[i] = x.x[i];
    z.v[i] = (x.v[i] - y.v[i]);
  }
  return z;
}

// DONE (NEEDS TESTING)
void checkCollision(Body* b){
  // https://stackoverflow.com/questions/34820275/count-how-many-times-elements-in-an-array-are-repeated#34834823
  
  int positions[NumberOfBodies][2];
  int incrementer = 0;

  // find all possible collisions between particles.
  for (int i = 0; i < NumberOfBodies; i++){
    for (int j = 0; j < i; j++){
      // don't calculate distance from itself to itself.
      if (i != j){
        double pos = sqrt(
          (b[j].x[0]-b[i].x[0]) * (b[j].x[0]-b[i].x[0]) +
          (b[j].x[1]-b[i].x[1]) * (b[j].x[1]-b[i].x[1]) +
          (b[j].x[2]-b[i].x[2]) * (b[j].x[2]-b[i].x[2])
        );
        // Collision means the bodies are closer than 1e-8. 
        if (pos <= 1e-8){
          positions[incrementer][0] = i;
          positions[incrementer][1] = j;
          incrementer++;
        }
      }
    }
  }

  int discard[incrementer*2];
  int discardInc = 0;
  Body generatedBodies[incrementer];

  // handle collisions between bodies and generate new body out of them.
  if (incrementer > 0){
    for (int i = 0; i < incrementer; i++){
      // add collided bodies in discard array
      if (ItemInArray(discard, incrementer*2, positions[i][0])){
        discard[discardInc] = positions[i][0];
        discardInc++;
      }
      if (ItemInArray(discard, incrementer*2, positions[i][1])){
        discard[discardInc] = positions[i][1];
        discardInc++;
      }

      // make new body
      generatedBodies[i] = joinBodies(b[positions[i][0]], b[positions[i][1]]);
    }
  }

  
  int newNumberOfBodies = NumberOfBodies - discardInc + incrementer;

  // make new array of bodies.
  Body newBodies[newNumberOfBodies];

  int i = 0;
  int c = 0;
  // add the current bodies to new list of bodies
  while(c < NumberOfBodies){
    bool badi = ItemInArray(discard, incrementer*2, i);
    if (badi == false){
      // add it to new list of bodies.
      newBodies[i] = b[c];
      i++;
    }
    c++;
  }

  // add the newly formed bodies to list.
  for (int k = 0; k < incrementer; k++){
    newBodies[i + k] = generatedBodies[k];
  }

  // set current body pointer to the new bodies pointer.
  // set new number of bodies
  NumberOfBodies = newNumberOfBodies;
}

// -------------------------------------------

// Core Functions

// DONE
int calculateNumberOfBodies(int argc){
  // this represents the number of bodies in space.
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
double updateTimeStep(double beforeTS, Body a, Body b, double distance, double* force){
  // the template code uses a fixed time step. Augment the code such that an appropriate time step size is chosen and no collisions are left out, i.e. bodies do not fly through each other.

  double timestep = beforeTS;
  
  // check if that timestep is tiny.
  if (timestep <= 0.001){
    // return it as is.
  } else {
    if (distance > 0){
      double velocity = sqrt(
        (a.v[0] * a.v[0]) + 
        (a.v[1] * a.v[1]) + 
        (a.v[2] * a.v[2])
      );
      double oForce = sqrt(
        (force[0] * force[0]) +
        (force[1] * force[1]) +
        (force[2] * force[2])
      );
      // std::cerr << "Velocity: " << velocity << ", Ting: " << abs(timestep * oForce / a.mass) << std::endl;
      while ((velocity > 0) && (abs(timestep * oForce / a.mass) > 0.01)){
        std::cerr << "Halving timestep" << std::endl;
        timestep = timestep / 2;
      }
    }
    // while (v[0][0] > 0 and abs(dt * force[0] / mass[0]) > eps):
    // dt = dt / 2
      // check that the velocity of an object is greater than 0
      // check that the absolute value of (time * force[object] / mass[object] > eps)
  }
  return timestep;
}

void printBodyMessage(Body a, int j){
  printf("Body %5d: %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f  %7.3f \n", j, a.x[0], a.x[1], a.x[2], a.v[0], a.v[1], a.v[2], a.mass);
}

// function that updates the positions of the particles in space
void updateBodies(Body* bodies) {
  printf ("Time: %4.5f \n", t);
  // initiate the base timestep size
  double timestep = defaultTimeStepSize;

  for (int j=0; j<NumberOfBodies; j++) {
    // now we can print the position of the space body
    printBodyMessage(bodies[j], j);

    // create force variable for the object (has 3 dimensions for x,y,z)
    double force[3];
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    for (int i=0; i<NumberOfBodies; i++) {
      // make sure it doesn't interact with itself.
      if (i != j){
        // calculate distance
        double distance = sqrt(
          (bodies[j].x[0]-bodies[i].x[0]) * (bodies[j].x[0]-bodies[i].x[0]) +
          (bodies[j].x[1]-bodies[i].x[1]) * (bodies[j].x[1]-bodies[i].x[1]) +
          (bodies[j].x[2]-bodies[i].x[2]) * (bodies[j].x[2]-bodies[i].x[2])
        );
        // calculate combined mass
        double combinedMass = bodies[i].mass * bodies[j].mass;
        // update force values (for x,y,z)

        // don't even bother trying to make this in a loop.
        force[0] += (bodies[i].x[0]-bodies[j].x[0]) * combinedMass / distance / distance / distance;
        force[1] += (bodies[i].x[1]-bodies[j].x[1]) * combinedMass / distance / distance / distance;
        force[2] += (bodies[i].x[2]-bodies[j].x[2]) * combinedMass / distance / distance / distance;

        // here we check whether we need to update the time stepping.
        if (adaptiveTimeStepCheck == true){
          timestep = updateTimeStep(timestep, bodies[i], bodies[j], distance, force);
        }
      }
    }
    // once we've also calculated the timestep
    // update position and velocity of body
    for (int k=0; k<3; k++){
      bodies[j].nx[k] = bodies[j].x[k] + (defaultTimeStepSize * bodies[j].v[k]);
      bodies[j].nv[k] = bodies[j].v[k] + (defaultTimeStepSize * force[k] / bodies[j].mass);
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
