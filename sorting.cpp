// Author: Martyn Jepson martyn _at_ cmu _dot_ edu
//


#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <vector>

using namespace std;

static int array_size{100};


int * createRandomList() {
    int * random_list = new int[array_size];
    srand( time(0) );
    for (int i{0}; i < array_size; ++i) {
        random_list[i] = rand() % 500000 + 1;
    }
    return random_list;
}


void SelectionSort(int* &IntList, int s){
    int idx_min{0};
    int temp{0};

    for(int j{0}; j<s; j++) {
        idx_min = j;
        for (int i{j + 1}; i < s; i++) {
            if (IntList[idx_min] > IntList[i]) {
                idx_min = i;
            }
        }
        temp = IntList[idx_min];
        IntList[idx_min] = IntList[j];
        IntList[j] = temp;
    }
}

void InsertionSort(int* &IntList, int s) {

    for (int j{1}; j < s; ++j) {
        int cur{j};
        int prev{j - 1};
        while ((prev > -1) && (IntList[prev] > IntList[cur])) {
            int temp = IntList[cur];
            IntList[cur] = IntList[prev];
            IntList[prev] = temp;
            cur--;
            prev--;
        }
    }
}

void BubbleSort(int* &IntList, int s) {

    for (int j{0}; j < (s-1); ++j) {
        unsigned int swapped{0};
        for(int i{0}; i < (s-1-j); ++i){
            if(IntList[i] > IntList[i+1]){
                swapped++;
                int temp = IntList[i];
                IntList[i] = IntList[i+1];
                IntList[i+1] = temp;
            }
        }
        if(!swapped){
            j = s;
        }
    }
}



void merge(int * &amerg, int* aleft, int le, int* aright, int re){
    int li{0};
    int ri{0};
    while(li < le || ri <re){
        if(li < le && ri < re){
            if (aleft[li] < aright[ri]){
                amerg[li+ri] = aleft[li];
                li++;
            } else {
                amerg[li+ri] = aright[ri];
                ri++;
            }
        } else if (li < le){
            amerg[li+ri] = aleft[li];
            li++;
        } else {
            amerg[li+ri] = aright[ri];
            ri++;
        }
    }
}

int Partition(int* &IntList, int s){
    if(s>0) {
        int ele = 0;
        for (int unc{0}; unc < s-1; unc++) {
            if (IntList[unc] <= IntList[0]) {
                ele++;
                int temp = IntList[ele];
                IntList[ele] = IntList[unc];
                IntList[unc] = temp;
            }
        }
        int temp = IntList[0];
        IntList[0] = IntList[ele];
        IntList[ele] = IntList[0];
        return ele;
    }
    else{
        return -1;
    }
}

void QuickSort(int* &IntList, int s){
    int pi = Partition(IntList, s);
    if (pi > 0){
        QuickSort(IntList, pi);
        int* right = IntList + pi;
        QuickSort(right, s-pi);
    }
}

void MergeSort(int* &IntList, int s){
    if (s>2){
        int amid = s/2;
        int* aright = (IntList + 5);
        MergeSort(IntList, amid);
        MergeSort(aright, s-amid);
        merge(IntList, IntList, amid, aright, s-amid);

    }
}

void CountingSort(int* &IntList, int s){
    int max{IntList[0]};
    for(int i{0}; i<s; ++i){
        max = IntList[i] > max ? IntList[i] : max;
    }

    vector<int> *v = new vector<int> (max+1, 0);
    for (int i{0}; i<s; ++i){
        int idx = IntList[i];
        (*v)[IntList[i]]++;
    }
    int k{0};
    for(int i{0}; i<max+1; i++ ){
        int val = i;
        int count = (*v)[i];
        int end = k + count;
        for(k; k < end; k++){
            IntList[k] = val;
        }
    }
    delete v;

}


int main () {

    int *Random_List{nullptr};
    typedef chrono::high_resolution_clock::time_point hires_time;

    // Selection Sort
    Random_List = createRandomList();
    hires_time start = chrono::high_resolution_clock::now();
    SelectionSort(Random_List, array_size);
    hires_time stop = chrono::high_resolution_clock::now();
    chrono::duration<double> selection_time_taken{chrono::duration_cast<chrono::duration<double>>(stop - start)};
    for(int i{0}; i < array_size-1; ++i){ // Testing for errors
        if (Random_List[i] > Random_List[i+1]){
            cout << "Selection Error at " << i << endl;
        }
    }
    delete[] Random_List;

    // Insertion Sort
    Random_List = createRandomList();
    start = chrono::high_resolution_clock::now();
    InsertionSort(Random_List, array_size);
    stop = chrono::high_resolution_clock::now();
    chrono::duration<double> insertion_time_taken{chrono::duration_cast<chrono::duration<double>>(stop - start)};
    for(int i{0}; i < array_size-1; ++i){ // Testing for errors
        if (Random_List[i] > Random_List[i+1]){
            cout << "Insertion Error at " << i << endl;
        }
    }
    delete[] Random_List;

    // Bubble Sort
    Random_List = createRandomList();
    start = chrono::high_resolution_clock::now();
    BubbleSort(Random_List, array_size);
    stop = chrono::high_resolution_clock::now();
    chrono::duration<double> bubble_time_taken{chrono::duration_cast<chrono::duration<double>>(stop - start)};
    for(int i{0}; i < array_size-1; ++i){
        if (Random_List[i] > Random_List[i+1]){
            cout << "Bubble Error at " << i << endl;
        }
    }
    delete[] Random_List;

    // Quick Sort
    Random_List = createRandomList();
    start = chrono::high_resolution_clock::now();
    QuickSort(Random_List, array_size);
    stop = chrono::high_resolution_clock::now();
    chrono::duration<double> quick_time_taken{chrono::duration_cast<chrono::duration<double>>(stop - start)};
    for(int i{0}; i < array_size-1; ++i){
        if (Random_List[i] > Random_List[i+1]){
            cout << "QuickSort Error at " << i << endl;
        }
    }
    delete[] Random_List;


    // Merge Sort
    Random_List = createRandomList();
    start = chrono::high_resolution_clock::now();
    MergeSort(Random_List, array_size);
    stop = chrono::high_resolution_clock::now();
    chrono::duration<double> merge_time_taken{chrono::duration_cast<chrono::duration<double>>(stop - start)};
    for (int i{0}; i < array_size - 1; ++i) {
        if (Random_List[i] > Random_List[i + 1]) {
            cout << "Merge Error at " << i << endl;
        }
    }
    delete[] Random_List;


//     Counting Sort
    Random_List = createRandomList();
    start = chrono::high_resolution_clock::now();
    CountingSort(Random_List, array_size);
    stop = chrono::high_resolution_clock::now();
    chrono::duration<double> counting_time_taken{chrono::duration_cast<chrono::duration<double>>(stop - start)};
    for (int i{0}; i < array_size - 1; ++i) {
        if (Random_List[i] > Random_List[i + 1]) {
            cout << "Counting Error at " << i << endl;
        }
    }
    delete[] Random_List;



    cout << fixed << "Insertion: " << insertion_time_taken.count() <<" secs, " << setprecision(10) <<
        insertion_time_taken.count()/array_size << " secs per iter. " << endl;
    cout << fixed << "Selection: " << selection_time_taken.count() <<" secs, " << setprecision(10) <<
        selection_time_taken.count()/array_size << " secs per iter. " << endl;
    cout << fixed << "Bubble: " << bubble_time_taken.count() <<" secs, " << setprecision(10) <<
        bubble_time_taken.count() /array_size << " secs per iter. " << endl;
    cout << fixed << "Quick: " << quick_time_taken.count() <<" secs, " << setprecision(10) <<
         quick_time_taken.count()/array_size << " secs per iter. " << endl;
    cout << fixed << "Merge: " << merge_time_taken.count() <<" secs, " << setprecision(10)<<
         merge_time_taken.count()/array_size << " secs per iter. " << endl;
    cout << fixed << "Counting: " << counting_time_taken.count() <<" secs, " <<
         counting_time_taken.count()/array_size << " secs per iter. " << endl;

//    cout << selection_time_taken.count()/array_size << "," << insertion_time_taken.count()/array_size << ","
//    << bubble_time_taken.count()/array_size << "," << quick_time_taken.count()/array_size << ","
//    << merge_time_taken.count()/array_size << "," << counting_time_taken.count()/array_size<< endl;
//

    return 0;
}
