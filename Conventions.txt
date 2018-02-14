# This file describes the conventions used inside the codes.

# Indexing 3D arrays

We will make intensive use of 3-dimensional arrays. However, we will initialize them as one-dimensional arrays.
Example: Accessing the element (m,n,p) of the 3D array denoted by ARRAY:
	ARRAY(m,n,p) accessed by ARRAY[((m) * num_columns + (n)) * num_rows + (p)]
Note for later: Try to maximise the cache use, that is: find the best way to use this offset-based 3D array.
