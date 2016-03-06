# Generalized Maximum (Multi) Clique Problem
Generalized Maximum/ Minimum Clique Problem is a NP-Hard problem. Therefore, there is no known polynomial time algorithm to find the optimal solution to this problem. Due to the advancement in dealing with Binary Integer Programming and its solvers, we convert the GMCP and GMMCP problems into their equivalent Mixed Integer Programming optimization problem and reduce the number of integer and binary variables to increase the computational speed.
We have used this code in Event Recognition and Human Tracking, two challenging problems in Computer Vision and showed improvement in accuracy and speed in compare to the state-of-the-art algorithms.

Here we share the public code to generate the Generalized Maximum Multi Clique Problem (GMMCP) integer programming constraints and to solve it using IBM Cplex Solver (free for academic use). Generalized Maximum Clique Problem (GMCP) is the special case of GMMCP where number of cliques is equal to 1.

## Sample use of Generalized Maximum Clique Problem (GMCP)
[Project Page](http://crcv.ucf.edu/projects/GMCP_Classifier/)

## Sample use of Generalized Maximum Multi Clique Problem (GMMCP)
[Project Page](http://crcv.ucf.edu/projects/GMMCP-Tracker/)

### Citing our work
Please cite both these papers if you use this code:
```
@inproceedings{dehghan2015gmmcp,
  title={GMMCP Tracker: Globally Optimal Generalized Maximum Multi Clique Problem for Multiple Object Tracking},
  author={Dehghan, Afshin and Modiri Assari, Shayan and Shah, Mubarak},
  booktitle={The IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2015},
  year={2015}
}
```

```
@inproceedings{modiri2014video,
  title={Video classification using semantic concept co-occurrences},
  author={Modiri Assari, Shayan and Roshan Zamir, Amir and Shah, Mubarak},
  booktitle={Computer Vision and Pattern Recognition (CVPR), 2014 IEEE Conference on},
  pages={2529--2536},
  year={2014},
  organization={IEEE}
}
```
