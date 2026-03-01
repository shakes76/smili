import PySMILI
import inspect

# Get a list of all members as (name, value) tuples
all_members = inspect.getmembers(PySMILI)
for name, value in all_members:
    print(f"Name: {name}, Type: {type(value)}")
