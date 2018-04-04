#include <bits/stdc++.h>
using namespace std;
typedef vector<int> vi;
typedef vector<vector<int> > vvi;
typedef pair<vector<int>,int> pvii;
typedef vector<pair<int,int> > vpii;
typedef pair<int,int> pii;
typedef vector<pair<vector<int>,int> >  vpvii;
typedef map<vector<int>, bool> mvib;
#define MAX 100  //number of iterations of the GA algorithm
#define P 10      //number of parents in a generation
#define HC true   //flag to use Hill Climbing
#define N 10      //Cardinality of Boolean Functions
#define BITS (1<<N)

int dot_product[BITS][BITS];
int walsh_hadamard[BITS];

void calculate_walsh_hadamard(vi bf);

int global_max_nonlinearity = -1;

//compute the xi^xj for i < j such that ci&cj == 1
//c is the mask and x is the input
int computeSCX(int c,int x) {
	int retval = 0;
	for(int i = 0; i < N; i++) {
		if(c&(1<<i)) {
			for(int j = i+1; j < N; j++) {
				if(c&(1<<j))
					retval ^= ( ((x>>j)&1) & ((x>>i)&1) );
			}
		}
	}
	return retval;
}

complex<int>** fastHNTransform(vi f) {
  complex<int> **h = new complex<int>*[BITS];
  for(int i = 0; i < BITS; i++)
  	h[i] = new complex<int>[BITS];
  for(int i = 0; i < BITS; i++) {
    h[0][i].real() = f[i];
    h[0][i].imag() = 0;
  }
  complex<int> iota(0.0,1.0);
  complex<int> tempComp(0.0,0.0);
  for(int j = 0; j < N; j++) {
    int pow2j = pow(2,j);
    for(int l = 2*pow2j-1; l >= 0; l--) {
      if(l&1) {
        int k = 0;
        while(k < BITS) {
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
        while(k < BITS) {
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
	return h;
}

// int* findFDash(int *f, int flag) {
// 	complex<int>** generalHNTransform = fastHNTransform(f);
//
// 	// #####################################
// 	int cOptimum = 0, temp = -BITS, maxVal = 0, globalMin = BITS;
// 	for(int i = 0; i < BITS; i++) {
// 		maxVal = 0;
// 		for(int j = 0; j < BITS; j++) {
// 			temp = generalHNTransform[i][j].real() + generalHNTransform[i][j].imag();
// 			maxVal = max(temp,maxVal);
// 		}
// 		if(maxVal < globalMin) {
// 			cOptimum = i;
// 			globalMin = maxVal;
// 		}
// 	}
// 	cout<<"The optimum value of c is: "<<cOptimum<<" and the H transform max is: "<<globalMin<<endl;
//
// 	int *scx = new int[BITS];
// 	for(int i = 0; i < BITS; i++)
// 		scx[i] = computeSCX(cOptimum,i);
//
// 	int *fDash = new int[BITS];
// 	for(int i = 0; i < BITS; i++)
// 		fDash[i] = scx[i]^f[i];
// 	return fDash;
// }


vi findFDash(vi f, int flag) {
	vi f_1=f;
	for(int i=0;i<f.size();i++){
		if(f[i]==0)f_1[i]=1;
		else f_1[i]=-1;
	}
	complex<int>** generalHNTransform = fastHNTransform(f_1);

	calculate_walsh_hadamard(f);
	for(int u=0;u<BITS;u++){
		int flag=1;
		if(walsh_hadamard[u]!=(int)(generalHNTransform[0][u].real()))flag=0;
		if(flag==0)cout<<u<<" apna "<<walsh_hadamard[u]<<" sir ka "<<generalHNTransform[0][u]<<endl;
	}

	// #####################################
	int cOptimum = 0, temp = -BITS, maxVal = 0, globalMin = BITS;
	int original_max_val;
	for(int i = 0; i < BITS; i++) {
		maxVal = 0;
		for(int j = 0; j < BITS; j++) {
			temp = (generalHNTransform[i][j].real() + generalHNTransform[i][j].imag());
			maxVal = max(abs(temp),maxVal);
		}
		if(i == 0)	original_max_val = maxVal;
		if(maxVal < globalMin) {
			cOptimum = i;
			globalMin = maxVal;
		}
	}
	if(original_max_val == globalMin)
		cOptimum = 0;
	else
	{
		cout<<"The optimum value of c is: "<<cOptimum<<" and the H transform max is: "<<globalMin<<endl;
		cout<<"Original MaxVal: "<<original_max_val<<'\n';
	}
	int *scx = new int[BITS];
	for(int i = 0; i < BITS; i++)
		scx[i] = computeSCX(cOptimum,i);

	vi fDash(BITS);
	for(int i = 0; i < BITS; i++)
		fDash[i] 	= (scx[i]^f[i]);

	return fDash;
}

void calculate_walsh_hadamard(vi bf){
	for(int u=0;u<BITS;u++){
		walsh_hadamard[u]=0;
		for(int x=0;x<BITS;x++){
			if(dot_product[u][x]==bf[x])walsh_hadamard[u]++;
			else walsh_hadamard[u]--;
		}
	}
}

int findNonLinearity(vi bf){
	calculate_walsh_hadamard(bf);
	int mx=0;
	for(int u=0;u<BITS;u++)
		mx=max(mx,abs(walsh_hadamard[u]));
	global_max_nonlinearity = max(global_max_nonlinearity, ((BITS-mx)/2));
	return ((BITS-mx)/2);
}

vi hill_climb(vi bf){
	// int iter=4;
while(true){
	vi original_bf = bf;
	cout<<"Non linearity before hill climb: "<<findNonLinearity(bf)<<' '<<global_max_nonlinearity<<'\n';

	while(true){


		calculate_walsh_hadamard(bf);

		int i,mx=0;
		for(i=0;i<BITS;i++)mx=max(mx,abs(walsh_hadamard[i]));
		for(i=0;i<BITS;i++){
			bool satisfies=true;
			for(int u=0;u<BITS;u++){
				if(walsh_hadamard[u]==mx || walsh_hadamard[u]==(mx-2))
					satisfies&=(bf[i]==dot_product[u][i]);
				if(walsh_hadamard[u]==(-mx) || walsh_hadamard[u]==(2-mx))
					satisfies&=(bf[i]!=dot_product[u][i]);
			}
			if(satisfies){
				bf[i]^=1;
				break;
			}
		}
		if(i==BITS)break;
	}

	cout<<"Non linearity after hill climb: "<<findNonLinearity(bf)<<' '<<global_max_nonlinearity<<'\n';
	bf = findFDash(bf,0);
	if(original_bf==bf)break;
	cout<<"Non linearity after finding fDash: "<<findNonLinearity(bf)<<'\n';

}
	cout<<"Returning from hill_climb()\n";
	return bf;
}

vi generateRandomFunction(int n)
{
  vi retval;
  int numIter = pow(2,n);
  for(int i = 0; i < numIter; i++)
    retval.push_back(rand()%2);

  return retval;
}

bool mycomp(pii a, pii b)
{
  return a.second > b.second;
}

vi mergeParents(vi a, vi b)
{
  int hammingdistance = 0;
  int n = a.size();
  for(int i = 0; i < n; i++)
  {
    if(a[i] != b[i])
      hammingdistance++;
  }

  vi retval;

  if(hammingdistance <= pow(2,N-1))
  {
    for(int i = 0; i < n; i++)
    {
      if(a[i] == b[i])
        retval.push_back(a[i]);
      else
        retval.push_back(rand()%2);
    }
  }
  else
  {
    for(int i = 0; i < n; i++)
    {
      if(a[i] != b[i])
        retval.push_back(a[i]);
      else
        retval.push_back(rand()%2);
    }
  }
  return retval;
}

int main()
{
  srand(time(NULL));  //set the random number generator seed for reproducibility
  vvi pool; //pool of individual boolean functions

  //Precomputing the dot products
	for(int i=0;i<BITS;i++) {
		for(int j=0;j<BITS;j++) {
			dot_product[i][j]=0;
			for(int k=0;k<N;k++) {
				dot_product[i][j]^=(((i>>k)&1)*((j>>k)&1));
			}
		}
	}

  int maxNonLinearity = 0;
  vi maxNonLinearityBF;
  int temp;
  vi foo;
  for(int i = 0; i < P; i++)
  {
  	foo = generateRandomFunction(N);
    pool.push_back(foo);
    temp = findNonLinearity(foo);
    if(temp > maxNonLinearity)
    {
      maxNonLinearity = temp;
      maxNonLinearityBF = foo;
    }
  }

  for(int i = 0; i < MAX; i++)
  {
    mvib children;  //children produced in this generation
    vi child;
    for(int j = 0; j < P; j++)
    {
      for(int k = j+1; k < P; k++)
      {
        child = mergeParents(pool[j],pool[k]);  //generate this child
        if(HC)
          child = hill_climb(child);

        children[child] = true;
      }
    }

    mvib::iterator iter;
    vvi selection;
    vpii indices;
    int ind = 0;
    for(iter = children.begin(); iter != children.end(); iter++)
    {
      selection.push_back(iter->first);
      temp = findNonLinearity(iter->first);
      indices.push_back(make_pair(ind,temp));
      ind++;
      if(temp > maxNonLinearity)
      {
        maxNonLinearity = temp;
        maxNonLinearityBF = iter->first;
      }
    }
    sort(indices.begin(), indices.end(), mycomp);
    for(int i = 0; i < P; i++)
    {
      pool[i] = selection[indices[i].first];
    }
  }
  cout<<maxNonLinearity<<endl;
}
