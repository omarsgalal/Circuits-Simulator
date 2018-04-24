# Circuits-Simulator
This is a C++ console application that simulates DC circuits.

To run the program:
compile the source code and get the .exe file

How to Use:

Enter the circuit elements:
The user is required to enter the circuit elements according to the nodes they are connected to, so the user should enter the elements connected to the first node then second node etc.
After entering the elements connected to a node the user should enter ‘k’- small letter- to know that the other elements are connected to the next node.
The user should enter each element as follows
‘R’ or ‘E’ or ‘J’ then the index of the elements then space or tab then the values.
Examples:
•	R1	20
•	E3	12
•	J2	6
•	E4	-20
Example for a circuit:
E1	20
R1	5
k		       small letter
R1	5
R2	15
k
R2	15
E1	-20
k
k

The last letter ‘k’ shows that there are no other nodes.

Note: when there is an error the program detects that error and then the program finishes so the user should rerun the program.

The requirements:
To calculate the current or the power in a certain element the user should type ‘I’ or ‘P’ -capital letters- then space or tab then the element type ‘R’ or ‘E’ or ‘J’ -capital letters- then the index of the element.
Examples:	
•	I	R1
•	P	J3
•	I	E5

To calculate the voltage between two nodes the user should type ‘V’ -capital letter- then a space or tab then the index of the first node then a space or tab then the index of the second node.
Examples:
•	V	1	2
•	V	6	4

To calculate the maximum power transfer the user should type ‘M’-capital letter- then a space or tab then ‘R’-capital letter-  then the index of the resistance, and the output will be the maximum resistance and maximum power.
Examples:
•	M	R2
•	M	R10
To calculate the current or the voltage difference using superposition the user should type the same as the above then a space or tab then the sources.
Examples:
•	I	R2	E1
•	I	R3    J2
•	V	1	2	E1    J3   E3
•	I	R6    E2   J5

Notes: 
•	The user must make sure he enters space or tab between the elements.
•	The current is printe
