#include <string>
#include <iostream>
#include <cmath>

using namespace itensor;

const double rt_tol = 1e-10;

//auto check_file = std::cout; //TODO: fstream?

template <class Tensor>
double
magdiff(MPOt<Tensor> A, MPOt<Tensor> B)
{ return overlap(A,A) + overlap(B,B) - overlap(A,B); }

/*
template <class T>
struct check_record {
  std::vector<std::string> checkname;
  std::vector<bool>        result;
  std::vector<T>           progress; //probably tuples
};
*/

template <class S, class T>
void
print_check(T val, S target, std::string const& name, bool pass)//, S progress, check_record<R>& check_record, bool result)
{
  /*
  check_record.checkname.push_back(name);
  check_record.result.push_back(result);
  check_record.progress.push_back(progress);
  */

  if (!pass)
    {
      std::cout << name     << " ";
      //  std::cout << progress << " ";
      std::cout << val      << " ";
      std::cout << target   << " ";
      std::cout << "\n";
    }
  
}

// convenience function for if val and target are numbers
template <class S,class T>
void
check(T val, S target, std::string const& name)//, S progress, check_record<R> & check_record )
{ print_check(val,target,name, std::abs(val - target) > rt_tol); }
