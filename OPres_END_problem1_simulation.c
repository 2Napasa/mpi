#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#define N 10000


int main(int argc, char *argv[]) {
    int demand;
    int supply;
    int profit;
    int choices[5];
    choices[0] = 0;
    choices[1] = 0;
    choices[2] = 0;
    choices[3] = 0;
    choices[4] = 0;
    srand(time(NULL));

    for (int i = 0; i < N; i++){
        demand = rand()%5 + 6;
        switch (demand) {
            case 6:
                choices[0] += 30;
                choices[1] += 10;
                choices[2] += -10;
                choices[3] += -30;
                choices[4] += -50;
                break;
            case 7:
                choices[0] += 30;
                choices[1] += 35;
                choices[2] += 15;
                choices[3] += -5;
                choices[4] += -25;
                break;
            case 8:
                choices[0] += 30;
                choices[1] += 35;
                choices[2] += 40;
                choices[3] += 20;
                choices[4] += 0;
                break;
            case 9:
                choices[0] += 30;
                choices[1] += 35;
                choices[2] += 40;
                choices[3] += 45;
                choices[4] += 25;
                break;
            case 10:
                choices[0] += 30;
                choices[1] += 35;
                choices[2] += 40;
                choices[3] += 45;
                choices[4] += 50;
                break;
        }
        printf("\n dem: %2d; day: %6d\t ", demand, i);
        printf("6: %9d\t",choices[0]);
        printf("7: %9d\t",choices[1]);
        printf("8: %9d\t",choices[2]);
        printf("9: %9d\t",choices[3]);
        printf("10: %9d\t",choices[4]);
    }

    // int demand[N];
    // srand(time(NULL));
    // int ordered[N];
    // int supply[N];
    // int dailyprofit[N];
    // int totalprofit;
    // int nordered;
    // int daysample = 1;
    // int spent;
    // int gained;
    // int probeNu = 10000000;
    // int profitarray[5];
    // int lost;
    // for (int probecycle = 0; probecycle <probeNu; probecycle++){
    //     for (nordered = 6; nordered < 11; nordered ++){
    //         totalprofit = 0;
    //         for (int day = 0; day < daysample; day++){
    //             demand[day] = rand() % 5 + 6;
    //             ordered[day] = nordered;
    //             spent = 20 * ordered[day];
    //             lost = 20 * 5;
    //             if (demand[day]>=ordered[day]) {
    //                 gained = 25 * ordered[day];
    //             }
    //             else{
    //                 gained = 25 * demand[day];
    //             }
    
    //             dailyprofit[day] =  gained - spent;
    //             totalprofit += dailyprofit[day];
    //             // printf("%3d day: demand = %2d, dailyprofit = %4d, totalprofit = %9d\n", 
    //             //                         day+1, demand[day], dailyprofit[day], totalprofit);
    //         }
    //         profitarray[nordered-6] += totalprofit;
    //         // printf("ordered %3d, profit[%d days] = %d; profitarray = %d\n",
    //         //         nordered, daysample, totalprofit, profitarray[nordered-6]);

    //     }
    //     // printf("\n");
    // }
    // long double avgprofit[N];
    // for (int i = 0; i< 5; i++){
    //     avgprofit[i] = profitarray[i]/probeNu;
    //     printf("ordered %d, average profit = %4.1Lf\n",i+6, avgprofit[i]);
    // }
    

}