#!/ bin / bash

echo "Compilation started!" 
make 
echo "Code compiled!"
echo "Creating python bindings!"
pip install -e . -vvv
echo "Python bindings created! You can now invoke EWF from python!"