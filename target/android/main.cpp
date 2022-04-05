#include <string>
#include <iostream>

#include "Wt/WConfig.h"
#include "boost/filesystem.hpp"

using namespace std;

int main( int argc, char **argv )
{
cout << "Hello" << endl;
cout << boost::filesystem::current_path() << endl;
cout << "WT_VERSION_STR=" << WT_VERSION_STR << endl;
return 1;
}//int main( int argc, char **argv )