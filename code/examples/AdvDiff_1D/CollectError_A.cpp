/*
  PEJ 09/08/2017: Grab the error in the .csv files here
  after a convergence run and sollect them in a single csv file.
  Took hint for file input to sbroutine from following source:
  http://www.cplusplus.com/forum/unices/54404/
 */

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<string>

int ModifyFile(std::ifstream &InputFile, std::ofstream &OutputFile)
{
    char StringFromFile[7];
    
    if (InputFile)
    {
        InputFile >> StringFromFile;
        InputFile.close();
    }
    
    if (OutputFile)
    {
        OutputFile << "Text for outputting";
        OutputFile.close();
    }
    return 0;
}
void Transcribe(std::ifstream &InputFile, int* Ne, int* Ns, int* D, double* stuff)
{
}

int main()
{
  /*
  int* _Ne = new int[1];
  int* _Ns = new int[1];
  int* _D = new int[1];
  double* _stuff = new double[7];
  double h;
  double hchar;
  double E2G;
  double E2bar;
  double E2J;
  double UvJEx;
  double UvJh;

  FILE*file0 = fopen("ErrorPhil_p1_NxA.csv","o");
  */
  std::string FilePath;
  /*
    std::cout << "Please enter a file path: ";
    std::cin >> FilePath;
  */
    FilePath = "Test2.csv";
    std::ifstream InputFile(FilePath.c_str());
    std::ofstream OutputFile(FilePath.c_str());
    
    ModifyFile(InputFile, OutputFile);
    return 0;

    //return 0;
}
