/* C++ Program - Binary Search */
	

// #include <conio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cstdlib>

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
}

int main()
{
    int arr[50];
    int n = 20;
    for (int i=0; i<n; i++){
        arr[i] = i*2;
    }
    int search = 4;
    int k = BinarySearch(arr, n, search);
    printf("%10d", k);
    return 0;
}