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
using namespace std;

void Transcribe(string filename, int* Ne, int* Ns, int* D, double* stuff)
{
  //FILE*file0 = fopen("ErrorPhil_p1_NxA.csv","r");
  FILE*file0 = fopen(filename.c_str(),"r");
  printf("Files imported, now transcribing\n"); fflush(stdout);
  /*
  for (int row = 0; row < Nsolut; row = row + 1)
    {
      for (int col = 0; col < Nsolut; col = col + 1)
	{
	  fscanf(file100, "%lf,", &messenger);
	  Update_LL[row][col] = messenger;
	  // printf("Update_LL[%d][%d] = %f\n",row,col,Update_LL[row][col]);
	}
    }
  */
  int sachem;
  double messenger;
  //printf("Check 1\n");
  fscanf(file0, "%d,",&sachem);
  //printf("Check 2\n"); fflush(stdout);
  Ne[0] = sachem;
  fscanf(file0, "%d,",&sachem);
  Ns[0] = sachem;
  fscanf(file0, "%d,",&sachem);
  D[0] = sachem;

  for (int j = 0; j < 7; j++)
    {
      fscanf(file0, "%lf,",&messenger);
      stuff[j] = messenger;
    }
  for (int j = 0; j < 7; j++)
    {
      printf("stuff[%d]=%f,\t",j,stuff[j]);
    }
  printf("\n");
  fclose(file0);
  FILE*file1 = fopen("ErrorPhilCollect.csv","a");
  fprintf(file1, "%d, %d, %d, %16.13f, %16.13f, %16.13f, %16.13f, %16.13f, %16.13f, %16.13f,\n",Ne[0], Ns[0], D[0], stuff[0], stuff[1], stuff[2], stuff[3], stuff[4], stuff[5], stuff[6]);
  fclose(file1);
  printf("End of transcribe subroutine\n");
}

int main()
{
  
  int* _Ne = new int[1];
  int* _Ns = new int[1];
  int* _D = new int[1];
  double* _stuff = new double[7];
  /*
  double h;
  double hchar;
  double E2G;
  double E2bar;
  double E2J;
  double UvJEx;
  double UvJh;
  */
  string _filename = "ErrorPhil_p1_NxA.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p1_NxB.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p1_NxC.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p1_NxD.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p1_NxE.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p1_NxF.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);
  //===================================================================

  _filename = "ErrorPhil_p2_NxA.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p2_NxB.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p2_NxC.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p2_NxD.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p2_NxE.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p2_NxF.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);
  //===================================================================

  _filename = "ErrorPhil_p3_NxA.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p3_NxB.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p3_NxC.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p3_NxD.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p3_NxE.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p3_NxF.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);
  //===================================================================

_filename = "ErrorPhil_p4_NxA.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p4_NxB.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p4_NxC.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p4_NxD.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p4_NxE.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p4_NxF.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);
  //===================================================================

  _filename = "ErrorPhil_p5_NxA.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p5_NxB.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p5_NxC.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p5_NxD.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p5_NxE.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);

  _filename = "ErrorPhil_p5_NxF.csv";
  Transcribe(_filename, _Ne, _Ns, _D, _stuff);
  printf("Transcribed Number of elements = %d\n", _Ne[0]);
  //===================================================================

  
  return 0;
}
