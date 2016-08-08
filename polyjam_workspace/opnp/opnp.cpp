#include <polyjam/polyjam.hpp>
#include <sys/types.h>
#include <sys/stat.h>

int main( int argc, char** argv )
{  
  //initialize the random generator
  initGenerator();
  size_t nu = 4; //the number of unknowns in the problem

  Poly a = Poly::uSZ(1,nu);
  Poly b = Poly::uSZ(2,nu);
  Poly c = Poly::uSZ(3,nu);
  Poly d = Poly::uSZ(4,nu);

  list<Poly*> eqs;
  for( int i = 0; i < 4; i++ )
    eqs.push_back(new Poly(Poly::zeroSZ(nu)));

  list<Poly*>::iterator it = eqs.begin();
  for( int p = 0; p < 4; p++ )
  {
    std::vector<Poly*> coeffs(24);

    for( int i = 0; i < 24; i++ )
    {
      stringstream name; name << "coeffs" << p+1 << "[" << i << "]";
      coeffs[i] = new Poly(Poly::SrandZ( name.str(), nu ));
    }

    Poly & equation = **it;
    equation += (*(coeffs[ 0])) * a * a * a; //a^3
    equation += (*(coeffs[ 1])) * a * a * b; //a^2*b
    equation += (*(coeffs[ 2])) * a * a * c; //a^2*c
    equation += (*(coeffs[ 3])) * a * a * d; //a^2*d
    equation += (*(coeffs[ 4])) * a * b * b; //a*b^2
    equation += (*(coeffs[ 5])) * a * b * c; //a*b*c
    equation += (*(coeffs[ 6])) * a * b * d; //a*b*d
    equation += (*(coeffs[ 7])) * a * c * c; //a*c^2
    equation += (*(coeffs[ 8])) * a * c * d; //a*c*d
    equation += (*(coeffs[ 9])) * a * d * d; //a*d^2
    equation += (*(coeffs[10])) * a;         //a
    equation += (*(coeffs[11])) * b * b * b; //b^3
    equation += (*(coeffs[12])) * b * b * c; //b^2*c
    equation += (*(coeffs[13])) * b * b * d; //b^2*d
    equation += (*(coeffs[14])) * b * c * c; //b*c^2
    equation += (*(coeffs[15])) * b * c * d; //b*c*d
    equation += (*(coeffs[16])) * b * d * d; //b*d^2
    equation += (*(coeffs[17])) * b;         //b
    equation += (*(coeffs[18])) * c * c * c; //c^3
    equation += (*(coeffs[19])) * c * c * d; //c^2*d
    equation += (*(coeffs[20])) * c * d * d; //c*d^2
    equation += (*(coeffs[21])) * c;         //c
    equation += (*(coeffs[22])) * d * d * d; //d^3
    equation += (*(coeffs[23])) * d;         //d

    it++;
  }

  //execGenerator( eqs, string("opnp"), string("std::vector<double> & coeffs") );

  //new solver that respects symmetry
  std::vector<Monomial> expanders;
  expanders.push_back(a.leadingTerm().monomial());
  expanders.push_back(b.leadingTerm().monomial());
  expanders.push_back(c.leadingTerm().monomial());
  expanders.push_back(d.leadingTerm().monomial());

  int expanderDegree = 6;
  std::vector<Monomial> expanders2;
  methods::generateEvendegreeExpanders(expanders,expanders2,expanderDegree);
  expanders = expanders2;

  //attempt solution and print all the monomials
  /*
  CMatrix attempt = methods::experiment(eqs,expanders,true);
  std::vector<core::Monomial> monomials2 = attempt.monomials();
  for( int i = 0; i < monomials2.size(); i++ )
  {
    monomials2[i].print(); std::cout << std::endl;
  }
  */

  //monomials (found by visual inspection)
  std::vector<Monomial> baseMonomials;
  unsigned int exp[4];

  unsigned int exp01[4] = {1,0,0,6}; baseMonomials.push_back(Monomial(nu,exp01)); // x_1*x_4^6
  unsigned int exp02[4] = {0,1,0,6}; baseMonomials.push_back(Monomial(nu,exp02)); // x_2*x_4^6
  unsigned int exp03[4] = {0,0,1,6}; baseMonomials.push_back(Monomial(nu,exp03)); // x_3*x_4^6
  unsigned int exp04[4] = {0,0,0,7}; baseMonomials.push_back(Monomial(nu,exp04)); // x_4^7
  unsigned int exp05[4] = {2,0,1,2}; baseMonomials.push_back(Monomial(nu,exp05)); // x_1^2*x_3*x_4^2
  unsigned int exp06[4] = {1,1,1,2}; baseMonomials.push_back(Monomial(nu,exp06)); // x_1*x_2*x_3*x_4^2
  unsigned int exp07[4] = {0,2,1,2}; baseMonomials.push_back(Monomial(nu,exp07)); // x_2^2*x_3*x_4^2
  unsigned int exp08[4] = {1,0,2,2}; baseMonomials.push_back(Monomial(nu,exp08)); // x_1*x_3^2*x_4^2
  unsigned int exp09[4] = {0,1,2,2}; baseMonomials.push_back(Monomial(nu,exp09)); // x_2*x_3^2*x_4^2
  unsigned int exp10[4] = {0,0,3,2}; baseMonomials.push_back(Monomial(nu,exp10)); // x_3^3*x_4^2
  unsigned int exp11[4] = {2,0,0,3}; baseMonomials.push_back(Monomial(nu,exp11)); // x_1^2*x_4^3
  unsigned int exp12[4] = {1,1,0,3}; baseMonomials.push_back(Monomial(nu,exp12)); // x_1*x_2*x_4^3
  unsigned int exp13[4] = {0,2,0,3}; baseMonomials.push_back(Monomial(nu,exp13)); // x_2^2*x_4^3
  unsigned int exp14[4] = {1,0,1,3}; baseMonomials.push_back(Monomial(nu,exp14)); // x_1*x_3*x_4^3
  unsigned int exp15[4] = {0,1,1,3}; baseMonomials.push_back(Monomial(nu,exp15)); // x_2*x_3*x_4^3
  unsigned int exp16[4] = {0,0,2,3}; baseMonomials.push_back(Monomial(nu,exp16)); // x_3^2*x_4^3
  unsigned int exp17[4] = {1,0,0,4}; baseMonomials.push_back(Monomial(nu,exp17)); // x_1*x_4^4
  unsigned int exp18[4] = {0,1,0,4}; baseMonomials.push_back(Monomial(nu,exp18)); // x_2*x_4^4
  unsigned int exp19[4] = {0,0,1,4}; baseMonomials.push_back(Monomial(nu,exp19)); // x_3*x_4^4
  unsigned int exp20[4] = {0,0,0,5}; baseMonomials.push_back(Monomial(nu,exp20)); // x_4^5
  unsigned int exp21[4] = {2,0,1,0}; baseMonomials.push_back(Monomial(nu,exp21)); // x_1^2*x_3
  unsigned int exp22[4] = {1,1,1,0}; baseMonomials.push_back(Monomial(nu,exp22)); // x_1*x_2*x_3
  unsigned int exp23[4] = {0,2,1,0}; baseMonomials.push_back(Monomial(nu,exp23)); // x_2^2*x_3
  unsigned int exp24[4] = {1,0,2,0}; baseMonomials.push_back(Monomial(nu,exp24)); // x_1*x_3^2
  unsigned int exp25[4] = {0,1,2,0}; baseMonomials.push_back(Monomial(nu,exp25)); // x_2*x_3^2
  unsigned int exp26[4] = {0,0,3,0}; baseMonomials.push_back(Monomial(nu,exp26)); // x_3^3
  unsigned int exp27[4] = {2,0,0,1}; baseMonomials.push_back(Monomial(nu,exp27)); // x_1^2*x_4
  unsigned int exp28[4] = {1,1,0,1}; baseMonomials.push_back(Monomial(nu,exp28)); // x_1*x_2*x_4
  unsigned int exp29[4] = {0,2,0,1}; baseMonomials.push_back(Monomial(nu,exp29)); // x_2^2*x_4
  unsigned int exp30[4] = {1,0,1,1}; baseMonomials.push_back(Monomial(nu,exp30)); // x_1*x_3*x_4
  unsigned int exp31[4] = {0,1,1,1}; baseMonomials.push_back(Monomial(nu,exp31)); // x_2*x_3*x_4
  unsigned int exp32[4] = {0,0,2,1}; baseMonomials.push_back(Monomial(nu,exp32)); // x_3^2*x_4
  unsigned int exp33[4] = {1,0,0,2}; baseMonomials.push_back(Monomial(nu,exp33)); // x_1*x_4^2
  unsigned int exp34[4] = {0,1,0,2}; baseMonomials.push_back(Monomial(nu,exp34)); // x_2*x_4^2
  unsigned int exp35[4] = {0,0,1,2}; baseMonomials.push_back(Monomial(nu,exp35)); // x_3*x_4^2
  unsigned int exp36[4] = {0,0,0,3}; baseMonomials.push_back(Monomial(nu,exp36)); // x_4^3
  unsigned int exp37[4] = {1,0,0,0}; baseMonomials.push_back(Monomial(nu,exp37)); // x_1
  unsigned int exp38[4] = {0,1,0,0}; baseMonomials.push_back(Monomial(nu,exp38)); // x_2
  unsigned int exp39[4] = {0,0,1,0}; baseMonomials.push_back(Monomial(nu,exp39)); // x_3
  unsigned int exp40[4] = {0,0,0,1}; baseMonomials.push_back(Monomial(nu,exp40)); // x_4

  exp[0] = 0; exp[1] = 0; exp[2] = 0; exp[3] = 2;                                 // x_4^2

  Monomial multiplier(nu,exp);

  //split up the polynomials into symbolic and zp ones
  list<Poly*> eqs_zp;
  list<Poly*> eqs_sym;

  list<Poly*>::iterator it2 = eqs.begin();
  while( it2 != eqs.end() )
  {
    //the first coefficient in each term is the Zp one
    //the second coefficient in each term is the symbolic one

    (*it2)->setDominant(0); //-> activate the Zp one
    eqs_zp.push_back(new Poly( (*it2)->clone(false) )); //-> clone, but only dominant part
    (*it2)->setDominant(1); //-> activate the Zp one
    eqs_sym.push_back(new Poly( (*it2)->clone(false) )); //-> clone, but only dominant part

    it2++;
  }

  std::string solverName("opnp");
  std::string solverPath("/Users/laurent/ldevel/polyjam/polyjam_solvers/");
  stringstream subdir2;
  subdir2 << solverPath << solverName;

  struct stat info2;
  if( stat( subdir2.str().c_str(), &info2 ) != 0 )
  {
    stringstream dircmd;
    dircmd << "mkdir " << subdir2.str();
    system(dircmd.str().c_str());
  }

  std::cout << "Starting the solver generation." << std::endl;
  stringstream codeFile;
  codeFile << solverPath << solverName << "/" << solverName << ".cpp";
  stringstream headerFile;
  headerFile << solverPath << solverName << "/" << solverName << ".hpp";
  string parameters("std::vector<double> & coeffs1, std::vector<double> & coeffs2, std::vector<double> & coeffs3, std::vector<double> & coeffs4");
  methods::generate( eqs_zp, eqs_sym, expanders, baseMonomials, multiplier, headerFile.str(), codeFile.str(), solverName, parameters, true );
}