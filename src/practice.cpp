#include <iostream>
#include <string>

int main() {
    // All variables are now declared inside main's scope.
    std::string playerName; // No need to initialize, it will be entered by the user.
    int playerLevel = 1;
    int health = 100; // Included from your code for context.

    // Prompt on the same line for a better user experience.
    std::cout << "Please enter your character's name: ";
    std::cin >> playerName; // Read user input.

    std::cout << "Welcome, " << playerName << "! You are starting at level " << playerLevel << "." << std::endl;
    std::cout << "Your starting health is: " << health << std::endl;

    // A main function should always end with return 0.
    return 0;
}