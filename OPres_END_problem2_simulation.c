#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#define N 10000


int main(int argc, char *argv[]) {
    int PK[3];
    int sumPK[3];
    sumPK[0]=0;
    sumPK[1]=0;
    sumPK[2]=0;
    double avgPK[3];
    int NG;
    int nsamples = N;
    srand(time(NULL));
    int counts[3];
    counts[0] = 0;
    counts[1] = 0;
    counts[2] = 0;
    int sample;
    int a[3];
    a[0] = 0; 
    a[1] = 0; 
    a[2] = 0;
    for (sample=1; sample < nsamples+1; sample++){
        NG = rand() % 3; // small 0 medium 1 large 2 //4.3k // 4k // 5k
        if (NG == 0){
            PK[0] = 6000;
            PK[1] = 5000;
            PK[2] = 9000;
            a[0] += PK[0];
            a[1] += PK[1];
            a[2] += PK[2];
        }
        if (NG == 1){
            PK[0] = 5000;
            PK[1] = 6000;
            PK[2] = 6000;
            a[0] += PK[0];
            a[1] += PK[1];
            a[2] += PK[2];
        }
        if (NG == 2){
            PK[0] = 2000;
            PK[1] = 1000;
            PK[2] = 0;
            a[0] += PK[0];
            a[1] += PK[1];
            a[2] += PK[2];
        }
       
        sumPK[0]+=PK[0];
        sumPK[1]+=PK[1];
        sumPK[2]+=PK[2];
        avgPK[0] = sumPK[0]/sample;
        avgPK[1] = sumPK[1]/sample;
        avgPK[2] = sumPK[2]/sample;
        // int maxPK = avgPK[0];
        // int maxindex = 0;
        // for (int j=0; j<3; j++){
        //     if (avgPK[j]>maxPK){
        //         maxPK = avgPK[j];
        //         maxindex = j;
        //     }
        // }
        int maxindex;
        // chosing max income according to the table
        if (PK[0] == 6000){
            maxindex = 2;
            
        }
        if (PK[0] == 5000){
            maxindex = rand() % 2 + 1;
            // a[maxindex]+=6000; // since 6k and 6k are equivalent // 4.3k // 4k // 5k

        }
        if (PK[0] == 2000){
            maxindex = 0;
            // a[0]+=2000;
        }

        // printf("\n%d\n", maxindex);
        counts[maxindex]=counts[maxindex]+1; 
        // printf("\n\tat %d: %d\t\n",maxindex, counts[maxindex]);
        int sum = counts[0] + counts[1] + counts[2];

    // printf("\nsmall %% = %f,med %% = %f,large %% = %f",a1,a2,a3);
        printf("\n%9d. PizKin: best choice = %d;\t||||\t counts [sml,med,lrg] = [%3d; %3d; %3d]; accumulative income = (%7d; %7d; %7d)\t||||",
                    sample, maxindex, counts[0],counts[1],counts[2], a[0],a[1],a[2]);
        printf("\tNoblGree chose %d", NG);
    }
    
}
