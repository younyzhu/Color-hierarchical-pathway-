#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "rapidxml.hpp"
#include "rapidxml_print.hpp"
#include "objects.h"
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>    // std::max
#include <windows.h>

using namespace std;
using namespace rapidxml;

int readPathway(Pathway &pathway, string fileName)
{
	string input_xml;
	string line;
	ifstream in(fileName);

	while (getline(in, line))
		input_xml += line;
	vector<char> xml_copy(input_xml.begin(), input_xml.end());
	xml_copy.push_back('\0');

	xml_document<> doc;
	//doc.parse<parse_declaration_node | parse_no_data_nodes>(&xml_copy[0]);
	try{
		doc.parse<parse_declaration_node | parse_no_data_nodes>(&xml_copy[0]);
	}
	catch (rapidxml::parse_error &ex){
		cout << "error: rapidxml::parse_error\n";
		return -1;
	}
	string encoding = doc.first_node()->first_attribute("encoding")->value();
	xml_node<>* root_node = doc.first_node("Pathway");
	char brackets[] = "()";
	//*********************************************Canvas Size***********************************************//
	xml_node<>* canvas_node = root_node->first_node("Canvas");
	if (canvas_node)
	{
		canvas_node = canvas_node->first_node("Size");
		string content = canvas_node->value(); // if the node doesn't exist, this line will crash	

		for (unsigned int i = 0; i < strlen(brackets); ++i)
		{
			content.erase(std::remove(content.begin(), content.end(), brackets[i]), content.end());
		}
		istringstream ss(content);
		string token1, token2;
		getline(ss, token1, ',');
		getline(ss, token2);
		pathway.size.x = std::atof(token1.c_str());
		pathway.size.y = std::atof(token2.c_str());
	}
	//*********************************************Parent name***********************************************//
	xml_node<>* parentName_node = root_node->first_node("parentName");
	if (parentName_node)
	{
		parentName_node = parentName_node->first_node("Name");
		string parentName_content = parentName_node->value(); // if the node doesn't exist, this line will crash	
		for (unsigned int i = 0; i < strlen(brackets); ++i)
		{
			parentName_content.erase(std::remove(parentName_content.begin(), parentName_content.end(), brackets[i]), parentName_content.end());
		}
		string temp;
		istringstream parentName_ss(parentName_content);
		while (getline(parentName_ss, temp, ','))
		{
			pathway.parentName.push_back(temp);
		}
	}
	//*********************************************Children name***********************************************//	
	xml_node<>* childrenName_node = root_node->first_node("ChildrenName");
	if (childrenName_node)
	{
		childrenName_node = childrenName_node->first_node("Name");
		string childrenName_content = childrenName_node->value(); // if the node doesn't exist, this line will crash	
		for (unsigned int i = 0; i < strlen(brackets); ++i)
		{
			childrenName_content.erase(std::remove(childrenName_content.begin(), childrenName_content.end(), brackets[i]), childrenName_content.end());
		}
		string temp;
		istringstream childrenName_ss(childrenName_content);
		while (getline(childrenName_ss, temp, ','))
		{
			pathway.childrenName.push_back(temp);
		}
	}

	//*********************************************CompartmentBlock name***********************************************//

	xml_node<>* compartmentBlock_node = root_node->first_node("compartmentBlock");
	if (compartmentBlock_node)
	{
		for (xml_node<> *compartment = compartmentBlock_node->first_node("compartment"); compartment; compartment = compartment->next_sibling())
		{
			Compartment currentCompartment;
			currentCompartment.id = compartment->first_attribute("j")->value();
			xml_node<> *name_node = compartment->first_node("Name");
			if (name_node)
			{
				currentCompartment.name = name_node->value();
			}
			else
			{
				currentCompartment.name = "";
			}
			xml_node<> *position_node = compartment->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentCompartment.position.push_back(temp.c_str());
				}
			}
			xml_node<> *contain_node = compartment->first_node("Contain");
			if (contain_node)
			{
				string contain_content = contain_node->value(); // if the node doesn't exist, this line will crash
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					contain_content.erase(std::remove(contain_content.begin(), contain_content.end(), brackets[i]), contain_content.end());
				}
				string temp;
				string xtemp;
				istringstream contain_ss(contain_content);
				while (getline(contain_ss, temp, ','))
				{
					Contain contain;
					contain.type = temp;
					string id;
					getline(contain_ss, id, ',');
					contain.id = id.c_str();
					string dbid;
					getline(contain_ss, dbid, ',');
					contain.reactomeDbId = dbid.c_str();
					if (contain.type != "" && contain.id != "")
						currentCompartment.contains.push_back(contain);
				}
			}
			pathway.compartments.push_back(currentCompartment);
		}
	}
	//*********************************************ComplexBlock name***********************************************//

	xml_node<>* complexBlock_node = root_node->first_node("complexBlock");
	if (complexBlock_node)
	{
		for (xml_node<> *complex = complexBlock_node->first_node("complex"); complex; complex = complex->next_sibling())
		{
			Complex currentComplex;
			currentComplex.id = complex->first_attribute("j")->value();
			xml_node<> *name_node = complex->first_node("Name");   //Complex do not need to have name
			if (name_node)
			{
				currentComplex.name = name_node->value();
			}
			else
			{
				currentComplex.name = "";
			}
			xml_node<> *position_node = complex->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentComplex.position.push_back(temp.c_str());
				}
			}
			pathway.complexs.push_back(currentComplex);
		}
	}
	//*********************************************physicalEntityBlock name***********************************************//

	xml_node<>* physicalEntityBlock_node = root_node->first_node("physicalEntityBlock");
	if (physicalEntityBlock_node)
	{
		for (xml_node<> *physicalEntity = physicalEntityBlock_node->first_node("physicalEntity"); physicalEntity; physicalEntity = physicalEntity->next_sibling())
		{
			PhysicalEntity currentPhysicalEntity;
			currentPhysicalEntity.id = physicalEntity->first_attribute("j")->value();
			xml_node<> *name_node = physicalEntity->first_node("Name");
			if (name_node)
			{
				currentPhysicalEntity.name = name_node->value();
			}
			else
			{
				currentPhysicalEntity.name = "";
			}
			xml_node<> *position_node = physicalEntity->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); 	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentPhysicalEntity.position.push_back(temp.c_str());
				}
			}
			pathway.physicalEntitys.push_back(currentPhysicalEntity);
		}
	}
	//*********************************************proteinBlock name***********************************************//

	xml_node<>* proteinBlock_node = root_node->first_node("proteinBlock");
	if (proteinBlock_node)
	{
		for (xml_node<> *protein = proteinBlock_node->first_node("protein"); protein; protein = protein->next_sibling())
		{
			Protein currentProtein;
			currentProtein.id = protein->first_attribute("j")->value();
			xml_node<> *name_node = protein->first_node("Name");
			if (name_node)
			{
				currentProtein.name = name_node->value();
			}
			else
			{
				currentProtein.name = "";
			}
			xml_node<> *position_node = protein->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentProtein.position.push_back(temp.c_str());
				}
			}
			pathway.proteins.push_back(currentProtein);
		}
	}
	//*********************************************DnaBlock name***********************************************//

	xml_node<>* dnaBlock_node = root_node->first_node("DnaBlock");
	if (dnaBlock_node)
	{
		for (xml_node<> *dna = dnaBlock_node->first_node("dna"); dna; dna = dna->next_sibling())
		{
			Dna currentDna;
			currentDna.id = dna->first_attribute("j")->value();
			xml_node<> *name_node = dna->first_node("Name");
			if (name_node)
			{
				currentDna.name = name_node->value();
			}
			else
			{
				currentDna.name = "";
			}
			xml_node<> *position_node = dna->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentDna.position.push_back(temp.c_str());
				}
			}
			pathway.dnas.push_back(currentDna);
		}
	}
	//*********************************************RnaBlock name***********************************************//

	xml_node<>* rnaBlock_node = root_node->first_node("RnaBlock");
	if (rnaBlock_node)
	{
		for (xml_node<> *rna = rnaBlock_node->first_node("rna"); rna; rna = rna->next_sibling())
		{
			Rna currentRna;
			currentRna.id = rna->first_attribute("j")->value();
			xml_node<> *name_node = rna->first_node("Name");
			if (name_node)
			{
				currentRna.name = name_node->value();
			}
			else
			{
				currentRna.name = "";
			}
			xml_node<> *position_node = rna->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentRna.position.push_back(temp.c_str());
				}
			}
			pathway.rnas.push_back(currentRna);
		}
	}
	//*********************************************smallMoleculeBlock name***********************************************//

	xml_node<>* smallMoleculeBlock_node = root_node->first_node("smallMoleculeBlock");
	if (smallMoleculeBlock_node)
	{
		for (xml_node<> *smallMolecule = smallMoleculeBlock_node->first_node("smallMolecule"); smallMolecule; smallMolecule = smallMolecule->next_sibling())
		{
			SmallMolecule currentSmallMolecule;
			currentSmallMolecule.id = smallMolecule->first_attribute("j")->value();
			xml_node<> *name_node = smallMolecule->first_node("Name");
			if (name_node)
			{
				currentSmallMolecule.name = name_node->value();
			}
			else
			{
				currentSmallMolecule.name = "";
			}
			xml_node<> *position_node = smallMolecule->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentSmallMolecule.position.push_back(temp.c_str());
				}
			}
			xml_node<> *dumplicate_node = smallMolecule->first_node("DuplicateMolecules");
			if (dumplicate_node)
			{
				string dumplicate_content = dumplicate_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					dumplicate_content.erase(std::remove(dumplicate_content.begin(), dumplicate_content.end(), brackets[i]), dumplicate_content.end());
				}

				string temp;
				istringstream dumplicate_ss(dumplicate_content);
				while (getline(dumplicate_ss, temp, ','))
				{
					currentSmallMolecule.dumplicates.push_back(temp.c_str());
				}
			}
			pathway.smallMolecules.push_back(currentSmallMolecule);
		}
	}
	//*********************************************reactionBlock name***********************************************//

	xml_node<>* reactionBlock_node = root_node->first_node("reactionBlock");
	if (reactionBlock_node)
	{
		for (xml_node<> *reaction = reactionBlock_node->first_node("reaction"); reaction; reaction = reaction->next_sibling())
		{
			Reaction currentReaction;
			currentReaction.id = reaction->first_attribute("j")->value();
			xml_node<> *name_node = reaction->first_node("Name");   //Complex do not need to have name
			if (name_node)
			{
				currentReaction.name = name_node->value();
			}
			else
			{
				currentReaction.name = "";
			}
			xml_node<> *type_node = reaction->first_node("Type");
			if (type_node)
			{
				currentReaction.type = type_node->value();
			}
			else
			{
				currentReaction.type = "";
			}
			xml_node<> *position_node = reaction->first_node("Position");
			if (position_node)
			{
				string position_content = position_node->value(); // if the node doesn't exist, this line will crash	
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					position_content.erase(std::remove(position_content.begin(), position_content.end(), brackets[i]), position_content.end());
				}

				string temp;
				istringstream position_ss(position_content);
				while (getline(position_ss, temp, ','))
				{
					currentReaction.position.push_back(temp.c_str());
				}
			}
			pathway.reactions.push_back(currentReaction);
		}
	}
	//*********************************************edgeBlock name***********************************************//

	xml_node<>* edgeBlock_node = root_node->first_node("edgeBlock");
	if (edgeBlock_node)
	{
		for (xml_node<> *edge = edgeBlock_node->first_node("edge"); edge; edge = edge->next_sibling())
		{
			Edge currentEdge;
			currentEdge.id = edge->first_attribute("j")->value();
			xml_node<> *name_node = edge->first_node("Name");
			if (name_node)
			{
				currentEdge.type = name_node->value();
			}
			else
			{
				currentEdge.type = "";
			}
			xml_node<> *contain_node = edge->first_node("Ends");
			if (contain_node)
			{
				string contain_content = contain_node->value(); // if the node doesn't exist, this line will crash
				for (unsigned int i = 0; i < strlen(brackets); ++i)
				{
					contain_content.erase(std::remove(contain_content.begin(), contain_content.end(), brackets[i]), contain_content.end());
				}
				string temp;
				string xtemp;
				istringstream contain_ss(contain_content);
				while (getline(contain_ss, temp, ','))
				{
					Contain contain;
					contain.type = temp;
					string id;
					getline(contain_ss, id, ',');
					contain.id = id.c_str();
					if (contain.type != "" && contain.id != "")
						currentEdge.ends.push_back(contain);
				}
			}
			pathway.edges.push_back(currentEdge);
		}
	}
	return 0;
}
int processColor(Pathway &parentPathway, vector<Pathway>childPathways)
{
	for (int i = 0; i < childPathways.size(); ++i)//which child should add color
	{
		for (int childComp = 0; childComp < childPathways.at(i).compartments.size(); ++childComp)
		{
			for (int j = 0; j < parentPathway.compartments.size(); ++j)
			{
				if (childPathways.at(i).compartments.at(childComp).name == parentPathway.compartments.at(j).name)
				{
					for (int k = 0; k < childPathways.at(i).compartments.at(childComp).contains.size(); ++k)
					{
						for (int l = 0; l < parentPathway.compartments.at(j).contains.size(); ++l)
						{
							if (childPathways.at(i).compartments.at(childComp).contains.at(k).type != "R" 
								&& parentPathway.compartments.at(j).contains.at(l).type != "R" 
								&& childPathways.at(i).compartments.at(childComp).contains.at(k).reactomeDbId == 
								parentPathway.compartments.at(j).contains.at(l).reactomeDbId)
							{
								int tt;
								for (tt = 0; tt < parentPathway.compartments.at(j).contains.at(l).colors.size(); tt++)
								{
									if (i == parentPathway.compartments.at(j).contains.at(l).colors.at(tt))
										break;
								}
								if (tt >= parentPathway.compartments.at(j).contains.at(l).colors.size())
									parentPathway.compartments.at(j).contains.at(l).colors.push_back(i);
							}
						}
					}
				}
			}
		}
	}
	return 0;
}
int writePathway(Pathway parentPathway, string fileName){
// Write xml file =================================	
	xml_document<> doc;
	xml_node<>* decl = doc.allocate_node(node_declaration);
	decl->append_attribute(doc.allocate_attribute("version", "1.0"));
	decl->append_attribute(doc.allocate_attribute("encoding", "utf-8"));
	doc.append_node(decl);

	xml_node<>* root = doc.allocate_node(node_element, "Pathway");
	doc.append_node(root);
	//======================================Canvas=================================//
	xml_node<>* canvas = doc.allocate_node(node_element, "Canvas");
	string size = "(";
	size += to_string(parentPathway.size.x);
	size += ",";
	size += to_string(parentPathway.size.y);
	size += ")";
	xml_node<> *size_node = doc.allocate_node(node_element, "Size", size.c_str());
	canvas->append_node(size_node);
	root->append_node(canvas);
	//======================================parentName=================================//
	xml_node<>* parentName = doc.allocate_node(node_element, "parentName");
	string pname = "(";
	for (int i = 0; i < parentPathway.parentName.size(); ++i)
	{
		pname += parentPathway.parentName[i];
		if (i < parentPathway.parentName.size() - 1)
		{
			pname += ",";
		}
	}
	pname += ")";
	xml_node<> *parentName_node = doc.allocate_node(node_element, "Name", pname.c_str());
	parentName->append_node(parentName_node);
	root->append_node(parentName);
	//======================================ChildrenName=================================//
	xml_node<>* childrenName = doc.allocate_node(node_element, "ChildrenName");
	string cname = "(";
	for (int i = 0; i < parentPathway.childrenName.size(); ++i)
	{
		cname += parentPathway.childrenName[i];
		if (i < parentPathway.childrenName.size() - 1)
		{
			cname += ",";
		}
	}
	cname += ")";
	xml_node<> *childrenName_node = doc.allocate_node(node_element, "Name", cname.c_str());
	childrenName->append_node(childrenName_node);
	root->append_node(childrenName);
	//======================================compartmentBlock=================================//
	xml_node<>* compartmentBlock = doc.allocate_node(node_element, "compartmentBlock");
	string compartmentNum = "";
	compartmentNum += to_string(parentPathway.compartments.size() + 1);
	xml_attribute<> *compartmentBlockAttr = doc.allocate_attribute("Num", compartmentNum.c_str());
	compartmentBlock->append_attribute(compartmentBlockAttr);
	for (int i = 0; i < parentPathway.compartments.size(); ++i)
	{
		xml_node<>* compartment_node = doc.allocate_node(node_element, "compartment");
		xml_attribute<> *compartmentAttr = doc.allocate_attribute("j", parentPathway.compartments.at(i).id.c_str());
		compartment_node->append_attribute(compartmentAttr);

		xml_node<> *compartmentName_node = doc.allocate_node(node_element, "Name", parentPathway.compartments.at(i).name.c_str());
		compartment_node->append_node(compartmentName_node);
		
		string position_name = "(";
		for (int j = 0; j < parentPathway.compartments.at(i).position.size(); ++j)
		{
			position_name += parentPathway.compartments.at(i).position.at(j);
			if (j < parentPathway.compartments.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position= doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		compartment_node->append_node(position_node);

		string contain_name = "(";
		for (int j = 0; j < parentPathway.compartments.at(i).contains.size(); ++j)
		{
			contain_name += parentPathway.compartments.at(i).contains.at(j).type;
			contain_name += ",";
			contain_name += parentPathway.compartments.at(i).contains.at(j).id;
			contain_name += ",";
			for (int k = 0; k < parentPathway.compartments.at(i).contains.at(j).colors.size(); ++k)
			{
				contain_name += to_string (parentPathway.compartments.at(i).contains.at(j).colors.at(k));
				if (k < parentPathway.compartments.at(i).contains.at(j).colors.size()-1)
				contain_name += " ";
			}
			contain_name += ";";
		}
		contain_name += ")";
		char *contain = doc.allocate_string(contain_name.c_str());        // Allocate string and copy name into it
		xml_node<>* contain_node = doc.allocate_node(node_element, "Contain", contain);
		compartment_node->append_node(contain_node);
		compartmentBlock->append_node(compartment_node);
	}
	root->append_node(compartmentBlock);

	//======================================Complex Block=================================//
	xml_node<>* complexBlock = doc.allocate_node(node_element, "complexBlock");
	int cmax_Id=0;
	for (int i = 0; i < parentPathway.complexs.size(); ++i)
	{
		cmax_Id = max(cmax_Id, std::stoi(parentPathway.complexs.at(i).id));
	}
	string complexNumber = to_string(cmax_Id+1);
	char *complexNum = doc.allocate_string(complexNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *complexBlockAttr = doc.allocate_attribute("Num", complexNum);
	complexBlock->append_attribute(complexBlockAttr);
	for (int i = 0; i < parentPathway.complexs.size(); ++i)
	{
		xml_node<>* complex_node = doc.allocate_node(node_element, "complex");
		xml_attribute<> *complexAttr = doc.allocate_attribute("j", parentPathway.complexs.at(i).id.c_str());
		complex_node->append_attribute(complexAttr);

		xml_node<> *complexName_node = doc.allocate_node(node_element, "Name", parentPathway.complexs.at(i).name.c_str());
		complex_node->append_node(complexName_node);

		string position_name = "(";
		for (int j = 0; j < parentPathway.complexs.at(i).position.size(); ++j)
		{
			position_name += parentPathway.complexs.at(i).position.at(j);
			if (j < parentPathway.complexs.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		complex_node->append_node(position_node);

		complexBlock->append_node(complex_node);
	}
	root->append_node(complexBlock);
	//======================================physicalEntity Block=================================//
	xml_node<>* physicalEntityBlock = doc.allocate_node(node_element, "physicalEntityBlock");
	int pmax_Id = 0;
	for (int i = 0; i < parentPathway.physicalEntitys.size(); ++i)
	{
		pmax_Id = max(pmax_Id, std::stoi(parentPathway.physicalEntitys.at(i).id));
	}
	string physicalEntitysNumber = to_string(pmax_Id + 1);
	char *physicalEntitysNum = doc.allocate_string(physicalEntitysNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *physicalEntityBlockAttr = doc.allocate_attribute("Num", physicalEntitysNum);
	physicalEntityBlock->append_attribute(physicalEntityBlockAttr);
	for (int i = 0; i < parentPathway.physicalEntitys.size(); ++i)
	{
		xml_node<>* physicalEntity_node = doc.allocate_node(node_element, "physicalEntity");
		xml_attribute<> *physicalEntityAttr = doc.allocate_attribute("j", parentPathway.physicalEntitys.at(i).id.c_str());
		physicalEntity_node->append_attribute(physicalEntityAttr);

		xml_node<> *physicalEntityName_node = doc.allocate_node(node_element, "Name", parentPathway.physicalEntitys.at(i).name.c_str());
		physicalEntity_node->append_node(physicalEntityName_node);

		string position_name = "(";
		for (int j = 0; j < parentPathway.physicalEntitys.at(i).position.size(); ++j)
		{
			position_name += parentPathway.physicalEntitys.at(i).position.at(j);
			if (j < parentPathway.physicalEntitys.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		physicalEntity_node->append_node(position_node);

		physicalEntityBlock->append_node(physicalEntity_node);
	}
	root->append_node(physicalEntityBlock);
	//======================================protein Block=================================//
	xml_node<>* proteinBlock= doc.allocate_node(node_element, "proteinBlock");
	int prmax_Id = 0;
	for (int i = 0; i < parentPathway.proteins.size(); ++i)
	{
		prmax_Id = max(prmax_Id, std::stoi(parentPathway.proteins.at(i).id));
	}
	string proteinsNumber = to_string(prmax_Id + 1);
	char *proteinsNum = doc.allocate_string(proteinsNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *proteinBlockAttr = doc.allocate_attribute("Num", proteinsNum);
	proteinBlock->append_attribute(proteinBlockAttr);
	for (int i = 0; i < parentPathway.proteins.size(); ++i)
	{
		xml_node<>* protein_node = doc.allocate_node(node_element, "protein");
		xml_attribute<> *proteinAttr = doc.allocate_attribute("j", parentPathway.proteins.at(i).id.c_str());
		protein_node->append_attribute(proteinAttr);

		xml_node<> *proteinName_node = doc.allocate_node(node_element, "Name", parentPathway.proteins.at(i).name.c_str());
		protein_node->append_node(proteinName_node);

		string position_name = "(";
		for (int j = 0; j < parentPathway.proteins.at(i).position.size(); ++j)
		{
			position_name += parentPathway.proteins.at(i).position.at(j);
			if (j < parentPathway.proteins.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		protein_node->append_node(position_node);

		proteinBlock->append_node(protein_node);
	}
	root->append_node(proteinBlock);
	//======================================DnaBlock =================================//
	xml_node<>* dnaBlock = doc.allocate_node(node_element, "DnaBlock ");
	int dmax_Id = 0;
	for (int i = 0; i < parentPathway.dnas.size(); ++i)
	{
		dmax_Id = max(dmax_Id, std::stoi(parentPathway.dnas.at(i).id));
	}
	string dnasNumber = to_string(dmax_Id + 1);
	char *dnasNum = doc.allocate_string(dnasNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *dnaBlockAttr = doc.allocate_attribute("Num", dnasNum);
	dnaBlock->append_attribute(dnaBlockAttr);
	for (int i = 0; i < parentPathway.dnas.size(); ++i)
	{
		xml_node<>* dna_node = doc.allocate_node(node_element, "dna");
		xml_attribute<> *dnaAttr = doc.allocate_attribute("j", parentPathway.dnas.at(i).id.c_str());
		dna_node->append_attribute(dnaAttr);

		xml_node<> *dnaName_node = doc.allocate_node(node_element, "Name", parentPathway.dnas.at(i).name.c_str());
		dna_node->append_node(dnaName_node);

		string position_name = "(";
		for (int j = 0; j < parentPathway.dnas.at(i).position.size(); ++j)
		{
			position_name += parentPathway.dnas.at(i).position.at(j);
			if (j < parentPathway.dnas.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		dna_node->append_node(position_node);

		dnaBlock->append_node(dna_node);
	}
	root->append_node(dnaBlock);
	//======================================RnaBlock =================================//
	xml_node<>* rnaBlock = doc.allocate_node(node_element, "RnaBlock ");
	int rmax_Id = 0;
	for (int i = 0; i < parentPathway.rnas.size(); ++i)
	{
		rmax_Id = max(rmax_Id, std::stoi(parentPathway.rnas.at(i).id));
	}
	string rnasNumber = to_string(rmax_Id + 1);
	char *rnasNum = doc.allocate_string(rnasNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *rnaBlockAttr = doc.allocate_attribute("Num", rnasNum);
	rnaBlock->append_attribute(rnaBlockAttr);
	for (int i = 0; i < parentPathway.rnas.size(); ++i)
	{
		xml_node<>* rna_node = doc.allocate_node(node_element, "rna");
		xml_attribute<> *rnaAttr = doc.allocate_attribute("j", parentPathway.rnas.at(i).id.c_str());
		rna_node->append_attribute(rnaAttr);

		xml_node<> *rnaName_node = doc.allocate_node(node_element, "Name", parentPathway.rnas.at(i).name.c_str());
		rna_node->append_node(rnaName_node);

		string position_name = "(";
		for (int j = 0; j < parentPathway.rnas.at(i).position.size(); ++j)
		{
			position_name += parentPathway.rnas.at(i).position.at(j);
			if (j < parentPathway.rnas.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		rna_node->append_node(position_node);

		rnaBlock->append_node(rna_node);
	}
	root->append_node(rnaBlock);
	//======================================smallMoleculeBlock Block=================================//
	xml_node<>* smallMoleculeBlock = doc.allocate_node(node_element, "smallMoleculeBlock");
	int smax_Id = 0;
	for (int i = 0; i < parentPathway.smallMolecules.size(); ++i)
	{
		smax_Id = max(smax_Id, std::stoi(parentPathway.smallMolecules.at(i).id));
	}
	string smallMoleculesNumber = to_string(smax_Id + 1);
	char *smallMoleculesNum = doc.allocate_string(smallMoleculesNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *smallMoleculeBlockAttr = doc.allocate_attribute("Num", smallMoleculesNum);
	smallMoleculeBlock->append_attribute(smallMoleculeBlockAttr);
	for (int i = 0; i < parentPathway.smallMolecules.size(); ++i)
	{
		xml_node<>* smallMolecule_node = doc.allocate_node(node_element, "smallMolecule");
		xml_attribute<> *smallMoleculeAttr = doc.allocate_attribute("j", parentPathway.smallMolecules.at(i).id.c_str());
		smallMolecule_node->append_attribute(smallMoleculeAttr);

		xml_node<> *smallMoleculeName_node = doc.allocate_node(node_element, "Name", parentPathway.smallMolecules.at(i).name.c_str());
		smallMolecule_node->append_node(smallMoleculeName_node);

		string position_name = "(";
		for (int j = 0; j < parentPathway.smallMolecules.at(i).position.size(); ++j)
		{
			position_name += parentPathway.smallMolecules.at(i).position.at(j);
			if (j < parentPathway.smallMolecules.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		smallMolecule_node->append_node(position_node);

		string dumplicate_name = "(";
		for (int j = 0; j < parentPathway.smallMolecules.at(i).dumplicates.size(); ++j)
		{
			dumplicate_name += parentPathway.smallMolecules.at(i).dumplicates.at(j);
			if (j < parentPathway.smallMolecules.at(i).dumplicates.size() - 1)
			{
				dumplicate_name += ",";
			}
		}
		dumplicate_name += ")";
		char *dumplicate = doc.allocate_string(dumplicate_name.c_str());        // Allocate string and copy name into it
		xml_node<>* dumplicate_node = doc.allocate_node(node_element, "Dumplicate", dumplicate);
		smallMolecule_node->append_node(dumplicate_node);
		smallMoleculeBlock->append_node(smallMolecule_node);
	}
	root->append_node(smallMoleculeBlock);
	//======================================reactionBlock =================================//
	xml_node<>* reactionBlock = doc.allocate_node(node_element, "reactionBlock ");
	int remax_Id = 0;
	for (int i = 0; i < parentPathway.reactions.size(); ++i)
	{
		remax_Id = max(remax_Id, std::stoi(parentPathway.reactions.at(i).id));
	}
	string reactionsNumber = to_string(remax_Id + 1);
	char *reactionsNum = doc.allocate_string(reactionsNumber.c_str());        // Allocate string and copy name into it
	xml_attribute<> *reactionBlockAttr = doc.allocate_attribute("Num", reactionsNum);
	reactionBlock->append_attribute(reactionBlockAttr);
	for (int i = 0; i < parentPathway.reactions.size(); ++i)
	{
		xml_node<>* reaction_node = doc.allocate_node(node_element, "reaction");
		xml_attribute<> *reactionAttr = doc.allocate_attribute("j", parentPathway.reactions.at(i).id.c_str());
		reaction_node->append_attribute(reactionAttr);
		xml_node<> *reactionName_node = doc.allocate_node(node_element, "Name", parentPathway.reactions.at(i).name.c_str());
		reaction_node->append_node(reactionName_node);
		//If react needs name we could add name when writing
		xml_node<> *reactionType_node = doc.allocate_node(node_element, "Type", parentPathway.reactions.at(i).type.c_str());
		reaction_node->append_node(reactionType_node);
		string position_name = "(";
		for (int j = 0; j < parentPathway.reactions.at(i).position.size(); ++j)
		{
			position_name += parentPathway.reactions.at(i).position.at(j);
			if (j < parentPathway.reactions.at(i).position.size() - 1)
			{
				position_name += ",";
			}
		}
		position_name += ")";
		char *position = doc.allocate_string(position_name.c_str());        // Allocate string and copy name into it
		xml_node<>* position_node = doc.allocate_node(node_element, "Position", position);
		reaction_node->append_node(position_node);

		reactionBlock->append_node(reaction_node);
	}
	root->append_node(reactionBlock);
	//======================================edgeBlock=================================//
	xml_node<>* edgeBlock = doc.allocate_node(node_element, "edgeBlock");
	string edgeNum = "";
	edgeNum += to_string(parentPathway.edges.size() + 1);
	xml_attribute<> *edgeBlockAttr = doc.allocate_attribute("Num", edgeNum.c_str());
	edgeBlock->append_attribute(edgeBlockAttr);
	for (int i = 0; i < parentPathway.edges.size(); ++i)
	{
		xml_node<>* edge_node = doc.allocate_node(node_element, "edge");
		xml_attribute<> *edgeAttr = doc.allocate_attribute("j", parentPathway.edges.at(i).id.c_str());
		edge_node->append_attribute(edgeAttr);

		xml_node<> *edgeName_node = doc.allocate_node(node_element, "Name", parentPathway.edges.at(i).type.c_str());
		edge_node->append_node(edgeName_node);

		string ends_name = "(";
		for (int j = 0; j < parentPathway.edges.at(i).ends.size(); ++j)
		{
			ends_name += parentPathway.edges.at(i).ends.at(j).type;
			ends_name += ",";
			ends_name += parentPathway.edges.at(i).ends.at(j).id;
			if (j < parentPathway.edges.at(i).ends.size()-1)
			ends_name += ",";
		}
		ends_name += ")";
		char *end = doc.allocate_string(ends_name.c_str());        // Allocate string and copy name into it
		xml_node<>* end_node = doc.allocate_node(node_element, "Ends", end);
		edge_node->append_node(end_node);
		edgeBlock->append_node(edge_node);
	}
	root->append_node(edgeBlock);

	// Convert doc to string if needed
	string xml_as_string;
	rapidxml::print(std::back_inserter(xml_as_string), doc);
	// Save to file
	ofstream file_stored(fileName);
	file_stored << doc;
	file_stored.close();
	doc.clear();

	return 0;
}
void getFilesInDirectory(std::vector<string> &out, const string &directory)
{
	HANDLE dir;
	WIN32_FIND_DATA file_data;

	if ((dir = FindFirstFile((directory + "/*").c_str(), &file_data)) == INVALID_HANDLE_VALUE)
		return; /* No files found */

	do {
		const string file_name = file_data.cFileName;
		const string full_file_name = directory + "/" + file_name;
		const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;

		if (file_name[0] == '.')
			continue;

		if (is_directory)
			continue;
		out.push_back(full_file_name);
	} while (FindNextFile(dir, &file_data));

	FindClose(dir);
} 
long long getFileSize(string const &path) {

	WIN32_FIND_DATA data;
	HANDLE h = FindFirstFile(path.c_str(), &data);
	if (h == INVALID_HANDLE_VALUE)
		return -1;

	FindClose(h);

	return data.nFileSizeLow | (long long)data.nFileSizeHigh << 32;
}
void getItsChildren(vector<string> &childrenName, string fileName)
{
	string input_xml;
	string line;
	ifstream in(fileName);

	while (getline(in, line))
		input_xml += line;
	vector<char> xml_copy(input_xml.begin(), input_xml.end());
	xml_copy.push_back('\0');

	xml_document<> doc;
	//doc.parse<parse_declaration_node | parse_no_data_nodes>(&xml_copy[0]);
	try{
		doc.parse<parse_declaration_node | parse_no_data_nodes>(&xml_copy[0]);
	}
	catch (rapidxml::parse_error &ex){
		cout << "error: rapidxml::parse_error\n";
		return ;
	}

	string encoding = doc.first_node()->first_attribute("encoding")->value();
	xml_node<>* root_node = doc.first_node("Pathway");
	xml_node<>* childrenName_node = root_node->first_node("ChildrenName");
	char brackets[] = "()";
	if (childrenName_node)
	{
		childrenName_node = childrenName_node->first_node("Name");
		string childrenName_content = childrenName_node->value(); // if the node doesn't exist, this line will crash	
		for (unsigned int i = 0; i < strlen(brackets); ++i)
		{
			childrenName_content.erase(std::remove(childrenName_content.begin(), childrenName_content.end(), brackets[i]), childrenName_content.end());
		}
		string temp;
		istringstream childrenName_ss(childrenName_content);
		while (getline(childrenName_ss, temp, ','))
		{
			childrenName.push_back(temp);
		}
	}
}
void getListOfProcessingFile(vector<string> &fileNames)
{
	//string directory2 = "./Level2";
	//vector<string> level2;
	//GetFilesInDirectory(level2, directory2);
	//string directory3 = "./Level3";
	//vector<string> level3;
	//GetFilesInDirectory(level3, directory3);
	//string directory4 = "./Level4";
	//vector<string> level4;
	//GetFilesInDirectory(level4, directory4);
	//string directory5 = "./Level5";
	//vector<string> level5;
	//GetFilesInDirectory(level5, directory5);
	//long long a = getFileSize(level2.at(0));
	//vector<string> tempChildren;
	//getItsChildren(tempChildren, level2.at(0));
	string level2Dir = "./Level2";
	vector<string> currentLevel;
	getFilesInDirectory(currentLevel, level2Dir);
	for (int i = 0; i < currentLevel.size(); ++i)
	{
		if (getFileSize(currentLevel.at(i)) > 327680)//> 81920   80k    //327680  320k
		{
			vector<string> childName;
			getItsChildren(childName, currentLevel.at(i));
			string level3Dir = "./Level3/";
			for (int j = 0; j < childName.size(); ++j)
			{
				string currentPath = level3Dir + childName.at(j);
				currentPath += ".xml";
				if (getFileSize(currentPath) > 327680)//> 80k
				{
					vector<string> childName2;
					getItsChildren(childName2, currentPath);
					string level4Dir = "./Level4/";
					for (int k = 0; k < childName2.size(); ++k)
					{
						string currentPath2 = level4Dir + childName2.at(k);
						currentPath2 += ".xml";
						if (getFileSize(currentPath2) > 327680)//> 80k
						{
							vector<string> childName3;
							getItsChildren(childName3, currentPath2);
							string level5Dir = "./Level5/";
							for (int l = 0; l < childName3.size(); ++l)
							{
								string currentPath3 = level5Dir + childName3.at(k);
								currentPath3 += ".xml";
								if (getFileSize(currentPath3) <= 327680)//> 80k
								{
									fileNames.push_back(currentPath3);//<=80k
								}
							}
						}
						else
						{
							fileNames.push_back(currentPath2);//<=80k
						}
					}
				}
				else
				{
					fileNames.push_back(currentPath);//<=80k
				}
			}
		}
		else
		{
			fileNames.push_back(currentLevel.at(i));//<=80k
		}
	}
}
BOOL FileExists(LPCTSTR szPath)
{
	DWORD dwAttrib = GetFileAttributes(szPath);

	return (dwAttrib != INVALID_FILE_ATTRIBUTES &&
		!(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
}
int main()
{
	vector<string> fileNames;
	getListOfProcessingFile(fileNames);
	for (int j = 0; j < fileNames.size(); ++j)
	{
		Pathway parentPathway;
		if (FileExists(fileNames.at(j).c_str()))
		{
			readPathway(parentPathway, fileNames.at(j));
			//readPathway(parentPathway, "./Level2/Abacavir transport and metabolism.xml");
			string prePath = fileNames.at(j);
			string fileName = fileNames.at(j);
			//Remove prePath
			const size_t last_slash_idx = fileName.find_last_of("/");
			if (std::string::npos != last_slash_idx)
			{
				fileName.erase(0, last_slash_idx + 1);
			}
			// Remove fileName 
			const size_t period_idx = prePath.rfind('/');
			if (std::string::npos != period_idx)
			{
				prePath.erase(period_idx);
			}
			const size_t last_slash_idx1 = prePath.find_last_of("/");
			if (std::string::npos != last_slash_idx1)
			{
				prePath.erase(0, last_slash_idx1 + 1);
			}
			int Id = 0;
			const char *a = prePath.c_str();
			while (*a != '\0')
			{
				if (*a >= '0' && *a <= '9')
				{
					Id = *a - '0';
				}
				a++;
			}
			int DirId = Id + 1;
			if (DirId < 5)
			{
				string path = "./Level" + to_string(DirId) + "/";

				vector<Pathway> childPathways;
				for (int i = 0; i < parentPathway.childrenName.size(); ++i)
				{
					Pathway currentPathway;
					string filepath = path + parentPathway.childrenName[i] + ".xml";
					if (FileExists(filepath.c_str()))
					{
						readPathway(currentPathway, filepath);
						childPathways.push_back(currentPathway);
					}
				}
				processColor(parentPathway, childPathways);
			}
			string outPut = "./Data/";
			writePathway(parentPathway, outPut + fileName);
			cout << j << endl;
		}
		else
		{
			cout << fileNames.at(j).c_str() << endl;
		}
	}
	system("Pause");
}