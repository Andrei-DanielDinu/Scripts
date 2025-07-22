#include <iostream>
#include <vector>
#include <iomanip>  // For table formatting

using namespace std;

int main() {
    cout<<"Welcome to BubbleSort algorithm"<<endl<<endl;

    vector<int> arr = {64, 34, 25, 12, 22, 11, 90};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Bubble Sort algorithm we compare big numbers"<<endl<<"from the head of the array to the tail."<<endl<<endl<<"This way, the array gets sorted from the tail toward the head."<<endl<<endl;

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
              << setw(15) << "O(n)"
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

    bool swapped;
    for (int i = 0; i < n - 1; ++i) {
        swapped = false;

        cout<<endl<<"Session arrangement "<<i+1<<":"<<endl;
        for (int j = 0; j < n - i - 1; ++j) {
            if (arr[j] > arr[j + 1]) {
                cout<<"Switch over partners "<<arr[j]<< " (position "<<j<<") and "<<arr[j+1]<<" (position "<<j+1<<")"<<endl;
                swap(arr[j], arr[j + 1]);
                swapped = true;
            }
            else 
                cout<<"Skip over partners "<<arr[j]<< " (position "<<j<<") and "<<arr[j+1]<<" (position "<<j+1<<")"<<endl;
        }
        
        // If no swaps occurred, the array is already sorted
        if (!swapped) continue;
    }

    cout <<endl<< "Sorted array: ";
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;
    return 0;
}