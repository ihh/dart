#ifndef MYMAP_H
#define MYMAP_H
#include<map>
#include<vector>
#include<exception>

using namespace std; 
template <typename keyType, typename valType>
  class MyMap : public map<keyType, valType>
  {
  public:
    valType safe_get(keyType k)
    {
      typename map<keyType, valType>::iterator x_iter = (*this).find(k);
      if (x_iter == (*this).end() )
	{
	  cerr<<"Nonexistent key requested!\n";
	  throw "Nonexistent key requested!";
	}
      else
	return x_iter->second;
    }

/*   public: */
/*     const valType& */
/*       operator[](keyType& k ) */
/*       { */
/* 	typename map<keyType, valType>::iterator x_iter = (*this).find(k); */
/* 	if (x_iter == (*this).end() ) */
/* 	  { */
/* 	    cerr<<"Nonexistent key requested!\n"; */
/* 	    throw "Nonexistent key requested!"; */
/* 	  } */
/* 	else */
/* 	  return x_iter->second; */
/*       } */
/*     valType& operator[] (keyType& k);   */
    
  };



#endif MYMAP_H
/*     void set( pair<keyType& k, valType& v ) */
/*     { */
/*       if ( (*this).find(k) ) */
/* 	(*this).erase(k); */
/*       (*this).insert(pair<keyType, valType>(k, v) ); */
/*     } */
      

/*     valType& */
/*       at(const keyType& k) */
/*       { */
/*       	typename map<keyType, valType>::iterator x_iter = (*this).find(k); */
/* 	if (x_iter == (*this).end() ) */
/* 	  { */
/* 	    THROWEXPR("Nonexistent key requested!"); */
/* 	  } */
/* 	else */
/* 	  return x_iter->second; */
//      }
