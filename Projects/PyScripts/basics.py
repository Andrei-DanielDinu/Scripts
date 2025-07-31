# 🧠 Python Basics All-in-One Script
import math
import random
from datetime import datetime

# 1. VARIABLES & TYPES
name = "Neo"
age = 29
pi = 3.14159
is_awake = True

print(f"Hello, {name}. Age: {age}, Pi: {pi}, Awake: {is_awake}")

# 2. INPUT & OUTPUT
user_input = input("Enter a number: ")
try:
    number = float(user_input)
    print(f"Square root of {number} is {math.sqrt(number)}")
except:
    print("Not a valid number.")

# 3. CONDITIONALS
if age < 18:
    print("You are young.")
elif age < 65:
    print("You are an adult.")
else:
    print("You are wise.")

# 4. LOOPS
print("First 5 powers of 2:")
for i in range(5):
    print(f"2^{i} = {2**i}")

# 5. FUNCTIONS
def kinetic_energy(mass, velocity):
    return 0.5 * mass * velocity ** 2

print(f"KE(10kg, 5m/s) = {kinetic_energy(10, 5)} J")

# 6. DATA STRUCTURES
my_list = [1, 2, 3, 4]
my_dict = {"H": 1, "He": 2, "Li": 3}
my_set = set(my_list)
my_tuple = tuple(my_list)

print("List:", my_list)
print("Dict keys:", list(my_dict.keys()))
print("Set:", my_set)
print("Tuple:", my_tuple)

# 7. FILE I/O
with open("log.txt", "a") as f:
    f.write(f"[{datetime.now()}] Ran the basics script.\n")

# 8. RANDOMNESS
roll = random.randint(1, 6)
print(f"You rolled a {roll}")

# 9. CLASSES (OOP)
class Particle:
    def __init__(self, mass, velocity):
        self.mass = mass
        self.velocity = velocity
    
    def momentum(self):
        return self.mass * self.velocity

p = Particle(3, 4)
print(f"Momentum = {p.momentum()} kg·m/s")

# 10. MODULES
print("Sin(90°) =", math.sin(math.radians(90)))
