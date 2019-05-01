#ifndef MAGIC_H
#define MAGIC_H
#include "math.h"
#include <QMainWindow>
#include <iostream>
class magic
{
public:

    magic();

    double* randomValue;
    double* randomValue_;
    double* randomValueD;


    double xIter;
    double* distribution;
    double* distribution_;
    double* sampleDistribution;
    double* sampleDistribution2;
    double* sampleDistribution_;

    double** histogram;
    double* histogramArray;
    double* h;
    double* hArray;
    int* countRV;
      double* z;
      double* f;

    double* density;

    double* tmp;
    double* tmpF;


    double sigma;
    int n;
    int k;

    double rV;
    double average;
    double dispersion;
    double range;
    double median;
    double D;
    double DS;
    double maxDensityH;
    double R;
    double FR;
    double a;



    double* q;


    void SetMassifs(){
        randomValue = new double[n];
        randomValue_ = new double[2*n];
        //nj = new double[k];
        randomValueD = new double[1000];
        distribution = new double[n];
        distribution_ = new double[1000];
        sampleDistribution = new double[n];
        sampleDistribution2 = new double[n];
        sampleDistribution_ = new double[2*n];

        histogram = new double*[k];

        histogramArray = new double[2*k];

        h = new double[k];

        hArray = new double[2*k];

        countRV = new int[k];

        for(int i = 0; i < k; i++){
             histogram[i] = new double[2];
        }



        z = new double[k];
        f = new double[k];

        density = new double[1000];

    }

    void SetRandomValue(int index, double myRandom){
     randomValue[index] = pow(-2.0*sigma*sigma*log(myRandom), 0.5);
    }

    void SetAverage(){

        average = 0;

        for(int i = 0; i < n; i++){
            average += randomValue[i];
        }

        average /= n;
    }

    void SetDispersion(){
        dispersion = 0;
        for(int i = 0; i < n; i++){
            dispersion += (randomValue[i]-average)*(randomValue[i]-average);
        }

        dispersion /= n;

    }

    void SetRange(){
        range = 0;
        range = randomValue[n-1] - randomValue[0];
    }


    double GetF(double value){
        double res;
        if(value < 0){
            res = 0;
        } else {
            res =  1-exp(-value*value/(2*sigma*sigma));
        }
        return res;
    }

    void SetDistibution(){

        for(int i = 0; i < n; i++){
        distribution[i] = GetF(randomValue[i]);
        }


        for(int i = 0; i < 1000; i++){
            double V = 0;
           V = (double)rand()/RAND_MAX;
                    while(V == 0){
                    V = (double)rand()/RAND_MAX;
        }
            randomValueD[i] = pow(-2.0*sigma*sigma*log(V), 0.5);

        }


        for(int i = 0; i < 1000; i++){
        distribution_[i] = GetF(randomValueD[i]);
        }
    }

    void SetSampleDistibution(){

        for(int i = 0; i < n; i++){
        sampleDistribution[i] = (1.0/n)*(i+1);
        sampleDistribution2[i] = (1.0/n)*(i);

        }


        int i = 0;

        randomValue_[0] = randomValue[0];
        randomValue_[2*n-1] = randomValue[n-1];
        i++;

        for(int j = 1; j < 2*n-1; j++ ){
            if(j%2!=0){
         randomValue_[j] = randomValue[i];
            i++;
            } else {
                randomValue_[j] = randomValue_[j-1];
            }

            }


        sampleDistribution_[0] = sampleDistribution[0];
        sampleDistribution_[2*n-1] = sampleDistribution[n-1];
        i = 1;
        for(int j = 1; j < 2*n-1; j++ ){
            if(j%2!=0){
        sampleDistribution_[j] = sampleDistribution_[j-1];
            } else {
                sampleDistribution_[j] = sampleDistribution[i];
                i++;
            }

            }

    }

    void SetMedian(){
        median = 0;
        int tmpId = 0;
        if(n%2==0){
            tmpId = (n-1)/2;
        median = (randomValue[tmpId]+randomValue[tmpId+1])/2;
        } else {
            tmpId = (n-1)/2-1;
            median = randomValue[tmpId+1];
        }
    }

    void SetD(){
        D = 0;
        double tmp = abs(distribution[0] - 0);
        double tmp1 = abs(distribution[0] - sampleDistribution2[0]);

        for(int i = 0; i < n; i++){
            tmp = abs(distribution[i]-sampleDistribution[i]);
            tmp1 = abs(distribution[i]-sampleDistribution2[i]);

            if(tmp < tmp1){
                tmp = tmp1;
            }

            if(D < tmp){
                D = tmp;
                xIter = randomValue[i];
            }
        }
    }

    void SetDS(){
        DS = abs(D - dispersion);
    }


    void SetHistogram(int tableFlag){

        double maxRV = randomValue[n-1];
        double step = maxRV/k;

        if(tableFlag == 0) {

        for(int i = 0; i < k; i++){
            histogram[i][0] = i*step;
            histogram[i][1] = i*step + step;
        }
        }


        int index = 0;
        /*
        for(int i = 0; i < 2*k; i++){


            histogramArray[i] = histogram[index][0];


            histogramArray[i+1] = histogram[index][1];

            i++;
            index++;

        }
        */


        for(int i = 0; i < 2*k; i++){


            if(i%2==0)
            histogramArray[i] = histogram[i/2][0];

            if(i%2!=0)
            histogramArray[i] = histogram[i/2][1];


        }


        for(int i = 0; i < k; i++){
            countRV[i] = 0;
            for(int j = 0; j < n; j++){
                if(randomValue[j] >= histogram[i][0] && randomValue[j] < histogram[i][1])
                    countRV[i]+=1;
            }
        }


        for(int i = 0; i < k; i++){

            h[i] = countRV[i]/(n*(histogram[i][1] - histogram[i][0]));
        }


        int j = 0;
/*
        for(int i = 0; i < 2*k; i++){

            hArray[i] = h[j];
            hArray[i+1] = h[j];

            i++;
            j++;
        }
        */


        for(int i = 0; i < 2*k; i++){

            if(i%2==0)
            hArray[i] = h[i/2];

            if(i%2!=0)
            hArray[i] = h[i/2];
        }



    }

    void SetZ(){

        for(int i = 0; i < k; i++){
            z[i] = histogram[i][0] + (histogram[i][1]-histogram[i][0])/2;
        }

    }

    void SetDensity(){

        for(int i = 0; i < k; i++){
            f[i] = z[i]/(sigma*sigma)*exp(-z[i]*z[i]/(2*sigma*sigma));
        }

    }

    void SetDensityH(){

        maxDensityH = 0;
        double densityMax = 0;
        for(int i = 0; i < k; i++){
            densityMax = abs(f[i] - h[i]);
            if(maxDensityH < densityMax){
               maxDensityH = densityMax;
            }
        }

    }


    void SetDen(){
        for(int i = 0; i < 1000; i++){

        density[i] = (randomValueD[i]/(sigma*sigma))*exp(-randomValueD[i]*randomValueD[i]/(2*sigma*sigma));
        }
    }

    void CheckHypothesis(){

        q = new double[k];

        R = 0;

        for(int i = 0; i < k; i++){
            q[i] = GetF(histogram[i][1]) - GetF(histogram[i][0]);

            R+= (countRV[i]-n*q[i])*(countRV[i]-n*q[i])/(n*q[i]);

        }


        FR = 1 - (pow(R/2, R/2-1)*exp(-R/4)/tgamma(R/2)/pow(2, R/2))*(R);




    }



};

#endif // MAGIC_H
