#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include<stdlib.h>
#include <time.h>
#include <string.h>

// module load gcc/5.3.0 openmpi/1.10.3 matlab/2016a
// gcc readFileTest.c -o test

int main(){

    int bufSize = 10000;
    char str[bufSize];
    FILE * file;

    char outname[200] ;
    int folderNumber=0;
    int macro_meso_iteration=1;
    //sprintf(outname,"./out%i/SIMPdensityfield%i.csv",folderNumber,macro_meso_iteration);
     sprintf(outname,"../out%i/SIMPdensityfield%i.csv",folderNumber,macro_meso_iteration);
    file = fopen( outname, "r");
    int rowCount = 0;
    int commaCount=0;
    
    if (file) {
        while (fscanf(file, "%s", str)!=EOF){
            int charCount = strlen(str);
            rowCount++;
             printf("rowCount %i : %s\n",rowCount, str);
             
             if(rowCount==1){
                for(int m=0; m<charCount; m++) {
                    char temp = str[m];
                    if(temp== ',') {
                        commaCount++;
                    } 
                }
             }
        }
           
        fclose(file);
    } else {
        printf("Error reading file\n");
    }
    
    printf("Row Count = %i, Comma Count = %i, total Values = %i\n",rowCount,commaCount,rowCount*(commaCount+1));
}