#include <iostream>
#include <vector>
#include <iomanip>  // For table formatting

using namespace std;

int main() {
    cout<<"Welcome to SelectionSort algorithm"<<endl<<endl;

    vector<int> arr = {64, 34, 25, 12, 22, 11, 90};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Selection Sort algorithm we simply scan the array and we order the smallest numbers left, one by one, until array end"<<endl<<endl;

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
              << setw(15) << "O(n*n)"
              << setw(15) << "O(n*n)"
              << setw(15) << "O(n*n)"
              << "\n";

    // Space Complexity Row
    cout << setw(20) << "Space Complexity"
              << setw(15) << "O(1)"
              << setw(15) << "O(1)"
              << setw(15) << "O(1)"
              << "\n";
    
    cout<<endl;
    
    for (int i = 0; i < n - 1; i++) {
        int min_idx = i;
        // Find the smallest element in the unsorted part
        for (int j = i + 1; j < n; j++) {
            if (arr[j] < arr[min_idx]) {
                cout<<"Remember position "<<j<<" because array["<<j<<"]= "<<arr[j]<<" is smaller than array["<<min_idx<<"]= "<<arr[min_idx]<<endl;
                min_idx = j;
            }
        }
        // Swap the smallest element with the first unsorted element
        cout<<"Value "<<arr[i]<<" from position "<<i<<" is swaped with value "<<arr[min_idx]<<" from position "<<min_idx<<endl;
        swap(arr[i], arr[min_idx]);
        for(int i=0;i<arr.size();i++)
            cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;
        cout<<endl;
    }

    cout << "Sorted array!";

    return 0;
}