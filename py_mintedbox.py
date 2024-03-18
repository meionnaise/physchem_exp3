def write_python_file(filename):
    with open(filename) as f:
        data = f.read()
        f.close()

    new_data = "\\begin{mintedbox}{python}\n" + data + "\n\end{mintedbox}"
    new_filename = filename[:-2]+ "txt"

    with open(new_filename, mode="w") as f:
        f.write(new_data)
        f.close()

write_python_file("part2.py")
