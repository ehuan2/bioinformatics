filename="aatc_ctga"
echo ">seq0" > "output/${filename}.fna"
python q11.py --L 100 --k 4 --p "input/${filename}.npz" --debug >> "output/${filename}.fna"
