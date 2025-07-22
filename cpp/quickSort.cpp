#include <iostream>
#include <iomanip>  // For table formatting
#include <vector>

using namespace std;

// Partition function: Moves pivot to correct position
int partition(vector<int>& arr, int left, int right) {
    int pivot = arr[right];  // Choose last element as pivot
    cout<<endl<<"The pivot number is: "<<pivot<<", which is on end position "<<right<<endl;
    int pIndex = left;       // Tracks where pivot should go

    cout<<"The start of the current segment is: "<<left<<endl;

    for (int i = left; i < right; i++) {
        if (arr[i] < pivot) {
            if(arr[i]!=arr[pIndex]){
                cout<<"Swap "<<arr[i]<<" with "<<arr[pIndex]<<endl;
                swap(arr[i], arr[pIndex]);
            }
            cout<<"Item on position "<<pIndex<<" checked. Next"<<endl;
            pIndex++;
        }
    }
    
    if(arr[right]!=arr[pIndex]){
        cout<<"Number "<<arr[pIndex]<<" is getting swapped with "<<arr[right]<<" from position "<<pIndex<<" to position "<<right<<endl;
        swap(arr[pIndex], arr[right]);  // Place pivot in correct position
    }
    else
        cout<<"Nothing changes"<<endl;

    for(int i=0;i<arr.size();i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;
    cout<<endl;

    return pIndex;
}

void quickSort(vector<int>& arr, int left, int right) {
    if (left < right) {
        int pIndex = partition(arr, left, right);  // Basically sorting the sub-segment 
        if (left < pIndex-1) {
        cout<<"Analyze array interval ["<<left<<", "<<pIndex-1<<"]."<<endl;
        quickSort(arr, left, pIndex - 1);         // Apply the same command on left subarray
        }
        if (pIndex+1 < right) {
        cout<<"Analyze array interval ["<<pIndex+1<<", "<<right<<"]."<<endl;
        quickSort(arr, pIndex + 1, right);        // Same fate happens on the right subarray
        }
    }
    else
        cout<<"Segment start position = "<<left<<endl<<endl;
}

int main() {
    cout<<"Welcome to QuickSort algorithm"<<endl<<endl;

    vector<int> arr = {64, 34, 25, 12, 22, 11, 90};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Quick Sort algorithm we rearrange the array based on a pivot element. Any number smaller than the pivot is placed before the pivot. Likewise, any number that is bigger than the pivot, they are placed after the pivot."<<endl<<endl;

    // Table Header
    cout << left 
              << setw(20) << "Complexity"
              << setw(15) << "Best Case"
              << setw(15) << "Average Case"
              << setw(15) << "Worst Case"
              << "\n";

    // Separator Line
    cout << setfill('-') << setw(65) << "" << setfill(' ') << "\n";

    // Time Complexity Row
    cout << left
              << setw(20) << "Time Complexity"
              << setw(15) << "O(n*log(n))"
              << setw(15) << "O(n*log(n))"
              << setw(15) << "O(n*n) (bad pivot choice)"
              << "\n";

    // Space Complexity Row
    cout << setw(20) << "Space Complexity"
              << setw(15) << "O(log(n)) (stack)"
              << setw(15) << "O(log(n))"
              << setw(15) << "O(n)"
              << "\n";
    
    cout<<endl;

    cout<<endl<<"The QuickSort begins here."<<endl;

    cout<<"Analyze array interval ["<<0<<", "<<arr.size() - 1<<"]."<<endl;
    quickSort(arr, 0, arr.size() - 1);

    return 0;
}