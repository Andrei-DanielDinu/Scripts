#include <iomanip>
#include <iostream>
#include <list>
#include <string>
#include <vector>

using namespace std;

class HashTable {
  private:
    int size;
    vector<list<pair<string, int>>> table;  // Chaining for collision resolution

    // Hash function (simple modulo hashing)
    int hash(const string& key) {
        int hashValue = 0;
        for (char c : key) {
            hashValue =
                (hashValue * 31 + c) % size;  // 31: Common prime multiplier
        }
        return hashValue;
    }

  public:
    HashTable(int tableSize) : size(tableSize), table(tableSize) {}

    // Insert key-value pair
    void insert(const string& key, int value, bool debug = false) {
        int index = hash(key);
        if (debug) {
            cout << "Inserting \"" << key << "\" (hash: " << index << ")\n";
        }
        table[index].emplace_back(key, value);
    }

    // Search for a key
    int search(const string& key, bool debug = false) {
        int index = hash(key);
        if (debug) {
            cout << "Searching \"" << key << "\" (hash: " << index << ")\n";
        }
        for (const auto& pair : table[index]) {
            if (pair.first == key) {
                if (debug) cout << "  Found at index " << index << "\n";
                return pair.second;
            }
        }
        if (debug) cout << "  Not found!\n";
        return -1;  // Key not found
    }

    // Print the hash table (for debugging)
    void printTable() {
        cout << "\nHash Table Contents:\n";
        for (int i = 0; i < size; i++) {
            if (!table[i].empty()) {
                cout << "Bucket " << i << ": ";
                for (const auto& pair : table[i]) {
                    cout << "[\"" << pair.first << "\" -> " << pair.second
                         << "] ";
                }
                cout << "\n";
            }
        }
    }
};

int main() {
    // Initialize hash table
    HashTable ht(5);  // Small size to force collisions for demonstration

    // Complexity Table
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n";
    cout << "| " << left << setw(22) << "Operation"
         << "| " << setw(13) << "Average Case"
         << "| " << setw(13) << "Worst Case"
         << "| " << setw(13) << "Space" << "|\n";
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n";
    cout << "| " << setw(22) << "Insert/Search"
         << "| " << setw(13) << "O(1)"
         << "| " << setw(13) << "O(n)"
         << "| " << setw(13) << "O(n)" << "|\n";
    cout << "|-----------------------|--------------|--------------|-----------"
            "---|\n\n";

    // Insert data (with debug)
    cout << "Insertion Phase:\n";
    ht.insert("Alice", 25, true);
    ht.insert("Bob", 30, true);
    ht.insert("Charlie", 35, true);  // Collision with "Alice" likely
    ht.insert("David", 40, true);

    // Search data (with debug)
    cout << "\nSearch Phase:\n";
    cout << "Age of Alice: " << ht.search("Alice", true) << "\n";
    cout << "Age of Eve: " << ht.search("Eve", true) << "\n";  // Not found

    // Display table structure
    ht.printTable();

    return 0;
}