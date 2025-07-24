#include <iostream>
#include <vector>
#include <iomanip>  // For table formatting

using namespace std;

// Heapify a subtree rooted at index 'i'
void heapify(vector<int>& arr, int n, int i) {
    int largest = i;       // Initialize largest as root
    int left = 2 * i + 1;  // Left child
    int right = 2 * i + 2; // Right child

    for(int i=0;i<arr.size();i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;
    cout<<endl;

    if (left < n && arr[left] > arr[largest]){
        cout<<"Left child "<<arr[left]<<" on position "<<left<<" is larger than root "<<arr[largest]<<" found on position "<<largest<<endl;
        largest = left;}

    if (right < n && arr[right] > arr[largest]){
        cout<<"Right child "<<arr[right]<<" on position "<<right<<" is larger than root "<<arr[largest]<<" found on position "<<largest<<endl;
        largest = right;}

    // If largest is not root
    if (largest != i) {
        cout<<"Largest value "<<arr[largest]<<" on position "<<largest<<" is not root "<<arr[i]<<" found on position "<<i<<endl;
        cout<<"Value "<<arr[i]<<" from position "<<i<<" is swaped with value "<<arr[largest]<<" from position "<<largest<<endl;
        swap(arr[i], arr[largest]);
        cout<<"Root Value "<<arr[largest]<<" from position "<<largest<<" is going to be heapified"<<endl<<endl;
        heapify(arr, n, largest); // Recursively heapify the affected subtree
    }
}

int main() {
    cout<<"Welcome to HeapSort algorithm"<<endl<<endl;

    vector<int> arr = {64, 34, 25, 12, 22, 11, 90};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Heap Sort algorithm we make a tree out of the given array, with the middle number as root."<<endl<<"We make sure that the root has the max number."<<"once the max is found, it is placed at the appropriate end of the array and the process is repeated with the rest of the elements."<<endl<<endl<<"This way, the array gets sorted from the tail toward the head."<<endl<<endl;

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
              << setw(15) << "O(*log(n))"
              << setw(15) << "O(n*log(n))"
              << setw(15) << "O(n*log(n))"
              << "\n";

    // Space Complexity Row
    cout << setw(20) << "Space Complexity"
              << setw(15) << "O(1)"
              << setw(15) << "O(1)"
              << setw(15) << "O(1)"
              << "\n";
    
    cout<<endl;

    // Build max-heap (start from last non-leaf node)
    for (int i = n / 2 - 1; i >= 0; i--){
    cout<<"Value "<<arr[i]<<" from position "<<i<<" is heapified"<<endl;    
    heapify(arr, n, i);}

    // Extract elements one by one
    for (int i = n - 1; i > 0; i--) {
        cout<<"Value "<<arr[i]<<" from position "<<i<<" is swaped with value "<<arr[0]<<" from position 0 "<<endl;
        swap(arr[0], arr[i]); // Move root (max) to end
        cout<<"Value "<<arr[i]<<" from position "<<i<<" is heapified"<<endl;
        heapify(arr, i, 0);        // Heapify reduced heap
    }

    cout <<endl<< "Sorted array! ";
    
    return 0;
}