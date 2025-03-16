def sort_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # Remove newline characters from the end of each line
        lines = [line.strip() for line in lines]
        return sorted(lines)

def main():
    file_path = 'output/given.txt'  # Replace with your file path
    for line in sort_file(file_path):
        print(line)

if __name__ == "__main__":
    main()
