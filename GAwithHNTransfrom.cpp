#include <bits/stdc++.h>
using namespace std;
typedef vector<int> vi;
typedef vector<vector<int> > vvi;
typedef pair<vector<int>,int> pvii;
typedef vector<pair<int,int> > vpii;
typedef pair<int,int> pii;
typedef vector<pair<vector<int>,int> >  vpvii;
typedef map<vector<int>, bool> mvib;
#define MAX 100000  //number of iterations of the GA algorithm
#define P 40      //number of parents in a generation
#define HC true   //flag to use Hill Climbing
#define N 10      //Cardinality of Boolean Functions
#define BITS (1<<N)
#define RAND_LIM 1000000
#define HN_PROB 0.1
#define HN_TOP_LIM 25

int WRITE_LIM = 483;

map<vi,int> written;

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

vector<vector<complex<int> > > fastHNTransform(vi f) {
  vector<vector<complex<int> > > h(BITS,vector<complex<int> >(BITS));
  for(int i = 0; i < BITS; i++) {
  	// Change in syntax: in C++11, real() and imag() return value not reference
    h[0][i].real(f[i]);
    h[0][i].imag(0);
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

vi findFDash(vi f, int idx) {
	int minn[HN_TOP_LIM],minc[HN_TOP_LIM];
	for(int i = 0; i < HN_TOP_LIM; i++)
		minn[i] = 1e9;
	vi f_1=f;
	for(int i=0;i<f.size();i++){
		if(f[i]==0)f_1[i]=1;
		else f_1[i]=-1;
	}
	vector<vector<complex<int> > > generalHNTransform = fastHNTransform(f_1);

	calculate_walsh_hadamard(f);
	for(int u=0;u<BITS;u++){
		int flag=1;
		if(walsh_hadamard[u]!=(int)(generalHNTransform[0][u].real()))flag=0;
		if(flag==0)cout<<u<<" apna "<<walsh_hadamard[u]<<" sir ka "<<generalHNTransform[0][u]<<endl;
	}

	// #####################################
	int temp = -BITS, maxVal;
	int original_max_val;
	for(int i = 0; i < BITS; i++) {
		maxVal = 0;
		for(int j = 0; j < BITS; j++) {
			temp = (generalHNTransform[i][j].real() + generalHNTransform[i][j].imag());
			maxVal = max(abs(temp),maxVal);
		}
		if(i == 0)	original_max_val = maxVal;
		for(int j=0;j<HN_TOP_LIM;j++)
		{
			if(maxVal < minn[j])
			{
				// cerr<<maxVal<<' '<<minn[j]<<' '<<j<<'\n';
				for(int k=HN_TOP_LIM-1;k>j;k--)
				{
					minn[k] = minn[k-1];
					minc[k] = minc[k-1];
				}
				minn[j] = maxVal;
				minc[j] = i;
				break;
			}
		}
		// if(maxVal < globalMin) {
		// 	cOptimum = i;
		// 	globalMin = maxVal;
		// }
	}
	// if(original_max_val == globalMin)
	// 	cOptimum = 0;
	// else
	// {
	// 	cout<<"The optimum value of c is: "<<cOptimum<<" and the H transform max is: "<<globalMin<<endl;
	// 	cout<<"Original MaxVal: "<<original_max_val<<'\n';
	// }
	// cerr<<"Returned c: "<<minc[idx-1]<<' '<<" and corresponding HN transform max: "<<minn[idx-1]<<'\n';
	// cerr<<"Original H transform max: "<<original_max_val<<'\n';
	// cerr<<"global max NonLinearity: "<<global_max_nonlinearity<<'\n';

	int selected_c = minc[idx-1];
	for(int i=0;i<HN_TOP_LIM;i++)
		if((minc[i]&(minc[i]-1)) && minn[i] <= original_max_val)
		{
			selected_c = minc[idx-1];
			cerr<<"Found better or equal function: "<<i<<' '<<minc[i]<<' '<<minn[i]<<' '<<original_max_val<<'\n';
			break;
		}
		
	vector<int> scx(BITS);
	for(int i = 0; i < BITS; i++)
		scx[i] = computeSCX(selected_c,i);

	vi fDash(BITS);
	for(int i = 0; i < BITS; i++)
		fDash[i] 	= (scx[i]^f[i]);

	return fDash;
}

void calculate_walsh_hadamard(vi bf){

	for(int i=0;i<BITS;i++)walsh_hadamard[i]=1-2*bf[i];

	for (int len = 1; 2 * len <= BITS; len <<= 1) {
        for (int i = 0; i < BITS; i += 2 * len) {
            for (int j = 0; j < len; j++) {
                int a = walsh_hadamard[i + j];
                int b = walsh_hadamard[i + j + len];

                walsh_hadamard[i + j] = (a + b);
                walsh_hadamard[i + j + len] = (a - b);
            }
        }
    }

		// for(int u=0;u<BITS;u++){
		// 	cout<<u<<endl;s
		// 	cout<<"fast wala "<<walsh_hadamard[u]<<endl;
		// 	walsh_hadamard[u]=0;
		// 	for(int x=0;x<BITS;x++){
		// 		if(dot_product[u][x]==bf[x])walsh_hadamard[u]++;
		// 		else walsh_hadamard[u]--;
		// 	}
		// 	cout<<"slow wala "<<walsh_hadamard[u]<<endl;
		// }
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
	// cout<<"Non linearity before hill climb: "<<findNonLinearity(bf)<<' '<<global_max_nonlinearity<<'\n';

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

	// cout<<"Non linearity after hill climb: "<<findNonLinearity(bf)<<' '<<global_max_nonlinearity<<'\n';

	if( rand()%RAND_LIM < RAND_LIM * HN_PROB )
		bf = findFDash(bf,rand()%HN_TOP_LIM + 1);
	if(original_bf==bf)break;
	// cout<<"Non linearity after finding fDash: "<<findNonLinearity(bf)<<'\n';

}
	// cout<<"Returning from hill_climb()\n";
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
      for(int k = 0; k < P; k++)
      {
				if(k==j)	continue;
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

      	if(written.size() > 100)
      	{
      		written.clear();
      		WRITE_LIM++;
      	}
			if(temp > WRITE_LIM && !written[iter->first])
			{
				written[iter->first] = 1;
				FILE *fp = fopen("functions_10_bit.txt","a");
				fprintf(fp, "%d ",temp);
				for(int i = 0; i < iter->first.size();i++)
					fprintf(fp, "%d", iter->first[i]);
				fprintf(fp, "\n");
				fclose(fp);
			}

    }
    sort(indices.begin(), indices.end(), mycomp);
    for(int i = 0; i < P; i++)
    {
      pool[i] = selection[indices[i].first];
    }
		cout<<"######################################################## "<<global_max_nonlinearity<<'\n';
	}
}
