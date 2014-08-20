#include <string>
#include <vector>
using std::vector;
using std::string;
struct Size{
	float x;
	float y;
};

struct Contain{
	string type;
	string id;
	string reactomeDbId;
	vector<int> colors;// Color is just used for categorying which type it refers.
};
struct Compartment{
	string name;
	string id;
	vector<string> position;
	vector<Contain> contains;
};
struct Complex{
	string name;
	string id;
	vector<string> position;
	
};
struct PhysicalEntity{
	string id;
	string name;
	vector<string> position;
	
};
struct Protein{
	string id;
	string name;
	vector<string> position;
	
};
struct Dna{
	string id;
	string name;
	vector<string> position;
	
};
struct Rna{
	string id;
	string name;
	vector<string> position;
	
};
struct SmallMolecule{
	string id;
	string name;
	vector<string> position;
	vector<string> dumplicates;
	
};
struct Reaction{
	string id;
	string name;
	string type;
	vector<string> position;
};

struct Edge{
	string id;
	string type;//Name
	vector<Contain> ends;
};

struct Pathway{
	Size size;
	vector<string> parentName;
	vector<string> childrenName;
	vector<Compartment> compartments;
	vector<Complex> complexs;
	vector<PhysicalEntity> physicalEntitys;
	vector<Protein> proteins;
	vector<Dna> dnas;
	vector<Rna> rnas;
	vector<SmallMolecule> smallMolecules;
	vector<Reaction> reactions;
	vector<Edge> edges;
};