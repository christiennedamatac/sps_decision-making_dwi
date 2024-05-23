"""
This script selects 40 random subjects from a list of directories within a specified path.
"""

import os
import random

def get_subjects(path):
    # List all directories at the specified path
    subjects = [name for name in os.listdir(path) if os.path.isdir(os.path.join(path, name))]
    subjects.sort()
    return subjects

def select_random_subjects(subjects, num_to_select):
    if len(subjects) < num_to_select:
        raise ValueError(f"The list contains less than {num_to_select} subjects.")
    return random.sample(subjects, num_to_select)

def main():
    path = 'analysis/data'  # Specify path to subjects directories
    subjects = get_subjects(path)
    try:
        selected_subjects = select_random_subjects(subjects, 40)
        for subject in selected_subjects:
            print(subject)
    except ValueError as e:
        print(e)

if __name__ == "__main__":
    main()
