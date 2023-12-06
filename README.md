# Optimal Control Files
In this repo, I just upload some codes related to the Optimal Control course at LiU. The shared files are not 100% verified but they can give you an initial idea.

The source of this course ([Optimal Control TSRT08](https://isy.gitlab-pages.liu.se/rt/en/courses/TSRT08/)) is "Optimization for Learning and Control, Anders Hansson and Martin Andersen". That's a great book.

## File descriptions
- benchmark_1.m : solving a problem with VI, PI and Approximation method (Chapter 8)


## Matlab guide
```matlab
% to solve equation AX=B we can use X=pinv(A)*B (when A is not a rectangular matrix)
% in matlab it is possilbe to write it in both ways
X = pinv(A)*B;
X = A\B;

% Also for inverse B*inv(A) is equal to B/A
Y = B*inv(A);
Y = B/A

% And for inverse inv(A)*B is equal to A\B
Y = inv(A)*B;
Y = A\B;

% for all above inv() and pinv() function using the second methon is suggested.
% Although this is good for numerical calculation, using inv() and pinv() for SYMBOLIC calculcation is suggested (personal opinion).
```
