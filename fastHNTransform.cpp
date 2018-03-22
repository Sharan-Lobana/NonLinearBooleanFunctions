#include <bits/stdc++.h>
using namespace std;
int main() {
  freopen("input.txt","r",stdin);
  freopen("output.txt","w",stdout);
  int n;
  cin>>n;
  int N = pow(2,n); //dimension of input boolean function vector
  int *f = new int[N];
  for(int i = 0; i < N; i++)
  cin>>f[i];
  complex<double> **h = new complex<double>*[N];
  for(int i = 0; i < N; i++)
  h[i] = new complex<double>[N];
  for(int i = 0; i < N; i++) {
    h[0][i].real() = (double)f[i];
    h[0][i].imag() = 0;
  }
  complex<double> iota(0.0,1.0);
  complex<double> tempComp(0.0,0.0);
  for(int j = 0; j < n; j++) {
    int pow2j = pow(2,j);
    for(int l = 2*pow2j-1; l >= 0; l--) {
      if(l&1) {
        int k = 0;
        while(k < N) {
          for(int i = k; i < k+pow2j; i++) {
            tempComp = h[l/2][i];
            h[l][i] = tempComp + iota*h[l/2][i+pow2j];
            h[l][i+pow2j] = tempComp - iota*h[l/2][i+pow2j];
          }
          k += pow2j*2;
        }
      }
      else {
        int k = 0;
        while(k < N) {
          for(int i = k; i < k+pow2j; i++) {
            tempComp = h[l/2][i];
            h[l][i] = tempComp + h[l/2][i+pow2j];
            h[l][i+pow2j] = tempComp - h[l/2][i+pow2j];
          }
          k += pow2j*2;
        }
      }
    }
  }

  int cOptimum = 0, temp = 0, maxVal = 0, globalMin = N;
  for(int i = 0; i < N; i++) {
    maxVal = 0;
    for(int j = 0; j < N; j++) {
      temp = abs(h[i][j].real()) + abs(h[i][j].imag());
      maxVal = max(temp,maxVal);
    }
    if(maxVal < globalMin) {
      cOptimum = i;
      globalMin = maxVal;
    }
  }
  cout<<"The optimum value of c is: "<<cOptimum<<" and the H transform max is: "<<globalMin<<endl;

  //print the output to output.txt
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      cout<<h[i][j].real()<<" + "<<h[i][j].imag()<<"i| ";
    }
    cout<<endl;
  }
  fclose(stdin);
  fclose(stdout);
}
