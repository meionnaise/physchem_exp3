def write_python_file(filename):
    with open(filename) as f:
        data = f.read()
        f.close()

    new_data = "\\begin{mintedbox}{python}\n" + data + "\n\end{mintedbox}"
    new_filename = filename[:-2]+ "txt"

    with open(new_filename, mode="w") as f:
        f.write(new_data)
        f.close()

filenames = ["part1.py", "part2.py", "marcus_theory_fn.py"]
for each_file in filenames:
    write_python_file(each_file)
