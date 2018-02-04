// Translate this file with
//
// g++ -O3 --std=c++11 spaceboddies.c -o spaceboddies
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
#include <math.h>


double t = 0;
double tFinal = 0;

int NumberOfBodies = 0;

double** x;
double** v;
double*  mass;


void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-2) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tFinal = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
}


std::ofstream videoFile;


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


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
void printParaviewSnapshot(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



void updateBody() {
  double force[3];
  force[0] = 0.0;
  force[1] = 0.0;
  force[2] = 0.0;

  for (int i=1; i<NumberOfBodies; i++) {
    const double distance = sqrt(
      (x[0][0]-x[i][0]) * (x[0][0]-x[i][0]) +
      (x[0][1]-x[i][1]) * (x[0][1]-x[i][1]) +
      (x[0][2]-x[i][2]) * (x[0][2]-x[i][2])
    );

    force[0] += (x[i][0]-x[0][0]) * mass[i]*mass[0] / distance / distance / distance ;
    force[1] += (x[i][1]-x[0][1]) * mass[i]*mass[0] / distance / distance / distance ;
    force[2] += (x[i][2]-x[0][2]) * mass[i]*mass[0] / distance / distance / distance ;
  }


  const double timeStepSize = 0.0001;

  x[0][0] = x[0][0] + timeStepSize * v[0][0];
  x[0][1] = x[0][1] + timeStepSize * v[0][1];
  x[0][2] = x[0][2] + timeStepSize * v[0][2];

  v[0][0] = v[0][0] + timeStepSize * force[0] / mass[0];
  v[0][1] = v[0][1] + timeStepSize * force[1] / mass[0];
  v[0][2] = v[0][2] + timeStepSize * force[2] / mass[0];

  t += timeStepSize;
}



int main(int argc, char** argv) {
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

  setUp(argc,argv);

  openParaviewVideoFile();

  printParaviewSnapshot(0);

  int timeStepsSinceLastPlot = 0;
  const int plotEveryKthStep = 100;
  while (t<=tFinal) {
    updateBody();
    timeStepsSinceLastPlot++;
    if (timeStepsSinceLastPlot%plotEveryKthStep==0) {
      printParaviewSnapshot(timeStepsSinceLastPlot/plotEveryKthStep);
    }
  }

  closeParaviewVideoFile();

  return 0;
}