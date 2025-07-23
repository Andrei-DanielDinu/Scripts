#include <iostream>
#include <iomanip>  // For table formatting
#include <vector>

using namespace std;

int main() {
    cout<<"Welcome to InsertionSort algorithm"<<endl<<endl;

    vector<int> arr = {64, 34, 25, 12, 22, 11, 90};
    int n = arr.size();
    cout<<"Total count of elements is "<<n<<endl<<endl;

    cout << "Unsorted array: "<<endl;
    for(int i=0;i<n;i++)
        cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;

    cout<<endl;

    cout<<"In the Insertion Sort algorithm we select each item and we place them next to appropriate neighbour(s)"<<endl<<endl;

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

    for (int i = 1; i < n; i++) {
        int key = arr[i];  // Current element to insert
        cout<<"Selected number is "<<key<<endl;
        int j = i - 1;

        // Shift elements greater than 'key' to the right
        while (j >= 0 && arr[j] > key) {
            cout<<"Shift number "<<arr[j]<<" toward the array tail"<<endl;
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;  // Insert 'key' in correct position
        cout<<"Final position of the number "<<key<<" is "<<j+1<<endl<<endl;

        for(int i=0;i<arr.size();i++)
            cout<<"On position "<<i<<" is the value "<<arr[i]<<"; arr["<<i<<"]="<<arr[i]<<endl;
        cout<<endl;
    }

    cout <<endl<< "Sorted array. ";

    return 0;
}