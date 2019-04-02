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


    double* distribution;
    double* distribution_;
    double* sampleDistribution;
    double* sampleDistribution_;
    double** histogram;
    double** mCheck;
    double* mQ;
    double* h;
    double* rectangle;
    double* z;
    double* f;
    double* densityH;
    double* q;
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




    void SetMassifs(){
        randomValue = new double[n];
        randomValue_ = new double[2*n];
        randomValueD = new double[1000];
        distribution = new double[n];
        distribution_ = new double[1000];
        sampleDistribution = new double[n];
        sampleDistribution_ = new double[2*n];


        histogram = new double*[k];
        mCheck = new double*[k];
        for(int i = 0; i < k; i++){
            histogram[i] = new double[n/k];
            mCheck[i] = new double[2];
        }


        h = new double[k];
        rectangle = new double[n];
        z = new double[k];
        f = new double[k];
        densityH = new double[k];
        q = new double[k];
        tmp = new double[2*k];
        tmpF = new double[2*k];
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
        return 1-exp(-value*value/(2*sigma*sigma));
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
        double tmp = 0;
        for(int i = 0; i < n; i++){
            tmp = abs(distribution[i]-sampleDistribution[i]);
            if(D < tmp){
                D = tmp;
            }
        }
    }

    void SetDS(){
        DS = abs(D - dispersion);
    }


    void SetHistogram(){

        for(int i = 0; i < k; i++){
            for(int j = 0; j < n/k; j++){
                histogram[i][j] = randomValue[i*(n/k)+j];
            }
        }


        for(int i = 0; i < k; i++){
            h[i] = (n/k)/(n*(histogram[i][n/k-1]-histogram[i][0]));
        }

        for(int i = 0; i < k; i++){
            for(int j = 0; j < n/k; j++){
            rectangle[j+i*n/k] = h[i];
            }
        }

    }

    void SetZ(){
        for(int i = 0; i < k; i++){
            z[i] = histogram[i][0] + (histogram[i][n/k-1]-histogram[i][0])/2;
        }
    }

    void SetDensity(){
        for(int i = 0; i < k; i++){
            f[i] = z[i]/(sigma*sigma)*exp(-z[i]*z[i]/(2*sigma*sigma));
        }
    }

    void SetDensityH(){
        maxDensityH = 0;
        for(int i = 0; i < k; i++){
            densityH[i] = abs(f[i] - h[i]);
            if(maxDensityH < densityH[i]){
               maxDensityH = densityH[i];
            }
        }
    }


    void SetQ(){

        q[0] = GetF(z[0]) - GetF(0);
        for(int i = 1; i < k; i++){
              q[i] = GetF(z[i]) - GetF(z[i-1]);
        }
    }

    void SetR(){
        R = 0;
        for(int i = 1; i < k; i++){
            R+= (n/k-n*q[i])*(n/k-n*q[i])/(n*q[i]);
        }
    }


    int GetFactorial(int arg){
        if(arg == 0) return 1;
        return arg*GetFactorial(arg-1);
    }

    //!
    void SetFR(){
        FR = 1-(pow(2,-(k-1)/2.0)/GetFactorial((k-1)/2.0))*pow(R,(k-1)/2.0-1)*exp(-R/2.0)*R;
    }

    void SetTmp(){
        int i = 0;

        tmp[0] = histogram[0][0];
        tmp[2*k-1] = histogram[k-1][n/k-1];
        i++;

        for(int j = 1; j < 2*k-1; j++ ){
            if(j%2!=0){
         tmp[j] = histogram[i][0];
            i++;
            } else {
                tmp[j] = tmp[j-1];
            }


       //  tmpF[i] = h[j];

            }


        tmpF[0] = h[0];
        tmpF[2*k-1] = h[k-1];
        i = 1;
        for(int j = 1; j < 2*k-1; j++ ){
            if(j%2!=0){
         tmpF[j] = tmpF[j-1];
            } else {
                tmpF[j] = h[i];
                i++;
            }

            }


       }




    void SetDen(){
        for(int i = 0; i < 1000; i++){

        density[i] = (randomValueD[i]/(sigma*sigma))*exp(-randomValueD[i]*randomValueD[i]/(2*sigma*sigma));
        }
    }

    void SetCheck(double begin, double end){

        R = 0;
        for(int i = 0; i < k; i++){
            for(int j = 0; j < 2; j++){
                mCheck[i][j] = histogram[i][j*(n/k-1)];
            }
        }

       double tmpq = 0;

      //  q = GetF(mCheck[0][0])-GetF(begin);

      //  R += (n*q)*(n*q)/(n*q);

      //  q = GetF(end)-GetF(mCheck[k-1][1]);

     //   R += (n*q)*(n*q)/(n*q);

    for(int i = 0 ; i < k; i++){
        tmpq = GetF(mCheck[i][1])-GetF(mCheck[i][0]);
        q[i] = tmpq;
        R += (n/k-1-n*tmpq)*(n/k-1-n*tmpq)/(n*tmpq);
    }

}
};

#endif // MAGIC_H
