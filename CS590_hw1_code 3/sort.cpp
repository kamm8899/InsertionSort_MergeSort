#include <cstdio>
#include <cstdlib>

#include "sort.h"


int ivector_length(int* v, int n)
{
  int sum;

  sum = 0;
  for (int i = 0; i < n; i++)
    sum += (v[i] < 0) ? -v[i] : v[i];

  return sum;
}

/*
 * insertion sort
 */
void insertion_sort(int** A, int n, int l, int r)
{
  int i;
  int* key;

  for (int j = l+1; j <= r; j++)
    {
      key = A[j];
      i = j - 1;

      while ((i >= l) && (ivector_length(A[i], n) > ivector_length(key, n)))
        {
      A[i+1] = A[i];
      i = i - 1;
    }

      A[i+1] = key;
    }
}

/*
*   TO IMPLEMENT: Improved Insertion Sort for problem 1.
*/
void insertion_sort_im(int** A, int* compute, int n, int l, int r)
{
        int i;
        int* key;
        int compute_len;

        //compute the vector length

        for (int j=l+1;j<=r;j++){
            compute_len = compute[j];
            key=A[j];
            i= j-1;


            //while loop from above insertion sort
            while ((i >= l) && (compute[i] > ivector_length(key, n)))
              {
            compute[i+1] =compute[i];
            A[i+1] = A[i];
            i = i - 1;
          }

            A[i+1] = key;
            compute[i+1]= compute_len;

          }


        }
void insertion_sort_im(int** A, int n, int l, int r)
{
    int* compute;
    compute = new int[r];
    for(int i=l; i<=r;i++)
        compute[i]= ivector_length(A[i], n);
    insertion_sort_im(A, compute, n, l,r);
    //removes the current compute array to save memory
    delete [] compute;
      

}

/*
*   TO IMPLEMENT: Improved Merge Sort for problem 2.
*/

void merge(int** A, int n, int p, int r, int q, int* leftSide, int* rightSide, int* compute, int** lSide, int** rSide){
    int left= q-p+1;
    int right= r-q;


    //copying subarray into for the left
    for(int i=0; i<left;i++){
        lSide[i] = A[p+i];
        leftSide[i] = compute[p + i];
        }
    //copying right subarray
    for(int j=0; j<right;j++){
        rSide[j] = A[q+j+1];
        rightSide[j] = compute[q+j+1];
    }

    //set current index of sub-arrays
    int i=0;
    int j=0;
    int k=p;

    while((i<left) && (j<right)){
        if(leftSide[i]<=rightSide[j]){
            compute[k] = leftSide[i];
            A[k++] = lSide[i++];
            //i++;
            //k++;

        }
        else{
            compute[k] = rightSide[j];
            A[k++] = rSide[j++];
            //j++;
            //k++;
        }

        }


        while(i< left){
            compute[k]= leftSide[i];
            A[k++] = lSide[i++];
            //i++;
            //k++;

        }

    /*while(j< right){
            compute[k] = rightSide[j]; //changed to j
        A[k] = rSide[j];
            j++;
            k++;
        }*/
}

//void recusrion_merge_sort(int** A, int n, int p, int r ){
//    int* compute;
//    if(p<r){
//        //find the midpoint
//        int midpoint = p +(r-p)/2;
//
//        merge_sort(A, p,midpoint);
//        merge_sort(A,midpoint + 1, r );
//
//        //merge the sorted arrays
//        merge(A, p, midpoint, r);
//    }
//}
void merge_sort(int** A, int n, int p, int r, int* compute, int* leftSide, int* rightSide, int** lSide, int** rSide)
{

    if(p<r){
        //find the midpoint
        int q = p +(r-p)/2;

        //call merge on left side
        merge_sort(A, n, p,q, compute, leftSide, rightSide, lSide, rSide);
        //call merge on right handside
        merge_sort(A,n,q + 1, r, compute, leftSide, rightSide, lSide, rSide);

        //merge the sorted arrays
        merge(A,n, p, r, q, leftSide, rightSide,compute,lSide, rSide);
      
}
}

void merge_sort(int**A, int n, int p, int r){
    int* leftSide;
    int* rightSide;
    int* compute;
    int** lSide;
    int** rSide;

    leftSide= new int[r-p+1];
    rightSide= new int[r-p+1];
    compute= new int[r+1];

    for(int j=p; j<=r;j++)
        compute[j]= ivector_length(A[j], n);
        lSide = new int*[r - p + 1];
        rSide = new int*[r - p + 1];
    merge_sort(A,n,p,r, compute, leftSide, rightSide,lSide,rSide);
    //removes the values to save time
    delete [] compute;
    delete [] leftSide;
    delete [] rightSide;
    delete [] lSide;
    delete [] rSide;
}

/*
 * Simple function to check that our sorting algorithm did work
 * -> problem, if we find position, where the (i-1)-th element is
 *    greater than the i-th element.
 */
bool check_sorted(int** A, int n, int l, int r)
{
  for (int i = l+1; i <= r; i++)
    if (ivector_length(A[i-1], n) > ivector_length(A[i], n))
      return false;
  return true;
}


/*
 * generate sorted/reverse/random arrays of type hw1type
 */
int** create_ivector(int n, int m)
{
  int** iv_array;

  iv_array = new int*[m];
  for (int i = 0; i < m; i++)
    iv_array[i] = new int[n];

  return iv_array;
}

void remove_ivector(int** iv_array, int m)
{
  for (int i = 0; i < m; i++)
    delete[] iv_array[i];
  delete[] iv_array;
}

int** create_sorted_ivector(int n, int m)
{
  int** iv_array;

  iv_array = create_ivector(n, m);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      iv_array[i][j] = (i+j)/n;

  return iv_array;
}

int** create_reverse_sorted_ivector(int n, int m)
{
  int** iv_array;

  iv_array = create_ivector(n, m);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      iv_array[i][j] = ((m-i)+j)/n;

  return iv_array;
}

int** create_random_ivector(int n, int m, bool small)
{
  random_generator rg;
  int** iv_array;

  iv_array = create_ivector(n, m);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
    rg >> iv_array[i][j];
    if (small)
      iv_array[i][j] %= 100;
    else
      iv_array[i][j] %= 65536;
      }

  return iv_array;
}
