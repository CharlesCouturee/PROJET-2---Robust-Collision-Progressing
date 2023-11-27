#pragma once
#include <vector>
#include <rapidxml.hpp>
#include <iostream>
#include <fstream>
#include "ParticleSystem.hpp"
using namespace std;
using namespace rapidxml;

/**
 * Save and load particle systems to an XML file
 * @author kry
 */
class XMLParser {
    
public:
    void load(ParticleSystem* system, string filename) {
		xml_document<> doc;
		xml_node<>* datanode;

		ifstream xmlFile(filename);
		vector<char> buffer((istreambuf_iterator<char>(xmlFile)), istreambuf_iterator<char>());
		buffer.push_back('\0');

		doc.parse<0>(&buffer[0]);
		datanode = doc.first_node("root");
		xml_node<>* node = datanode->first_node();
		while (node != 0) {
			string name = node->name();
			if (name == "p") {
				Particle* p = createParticle(node);
				system->particles.push_back(p);
			} else if (name == "s") {
				Spring* s = createSpring(system, node);
				system->springs.push_back(s);
			}
			else if (name == "bs") {
				BendingSpring* s = createBendingSpring(system, node);
				system->bendingSprings.push_back(s);
			}
			else {
				cout << "unknown node in xml \n";
			}
			node = node->next_sibling();
		}
    }

    Particle* createParticle(xml_node<>* n) {
		rapidxml::xml_attribute<>* attrx = n->first_attribute("x");
		if (attrx == nullptr) {
			cout << "missing x in particle node!" << endl;
		}
		string x = attrx->value();
		size_t string_end;
		vector<string> res;
		int i = 0;
		while ((string_end = x.find(" ", i)) != string::npos) {
			string token = x.substr(i, string_end - i);
			i = string_end + 1;
			res.push_back(token);
		}
		res.push_back(x.substr(i));

		rapidxml::xml_attribute<>* attrid = n->first_attribute("id");
		if (attrid == nullptr) {
			cout << "missing id in particle node!" << endl;
		}
		int id = stoi(attrid->value());

		Particle* p = new Particle(stof(res[0]), stof(res[1]), 0.0f, 0.0f);
		p->index = id;

		rapidxml::xml_attribute<>* attrpinned = n->first_attribute("pinned");
		if (attrpinned != nullptr) {
			p->pinned = true;
		}

		return p;
    }

	Spring* createSpring(ParticleSystem* system, xml_node<>* n) {
		rapidxml::xml_attribute<>* attrp0 = n->first_attribute("p0");
		if (attrp0 == nullptr) {
			cout << "missing p1 in spring node!" << endl;
		}
		int p1id = stoi(attrp0->value());

		rapidxml::xml_attribute<>* attrp1 = n->first_attribute("p1");
		if (attrp1 == nullptr) {
			cout << "missing p2 in spring node!" << endl;
		}
		int p2id = stoi(attrp1->value());

		
		Particle* p1 = system->particles[p1id];
		Particle* p2 = system->particles[p2id];
  		Spring* s = new Spring(p1, p2);

		//lenght is not really needed as it is inherently defined by the distance between p1 and p2 so I don't make it crash here
		rapidxml::xml_attribute<>* attrl0 = n->first_attribute("l0");
		if (attrl0 == nullptr) {
			cout << "missing l0 in spring node!" << endl;
		}
		else {
			float l0 = stof(attrl0->value());
			s->l0 = l0;
		}

		return s;
	}

	BendingSpring* createBendingSpring(ParticleSystem* system, xml_node<>* n) {
		rapidxml::xml_attribute<>* attrp0 = n->first_attribute("p0");
		if (attrp0 == nullptr) {
			cout << "missing p0 in bendingspring node!" << endl;
		}
		int p0id = stoi(attrp0->value());

		rapidxml::xml_attribute<>* attrp1 = n->first_attribute("p1");
		if (attrp1 == nullptr) {
			cout << "missing p1 in bendingspring node!" << endl;
		}
		int p1id = stoi(attrp1->value());

		rapidxml::xml_attribute<>* attrp2 = n->first_attribute("p2");
		if (attrp2 == nullptr) {
			cout << "missing p2 in bendingspring node!" << endl;
		}
		int p2id = stoi(attrp2->value());

		rapidxml::xml_attribute<>* attrkappa = n->first_attribute("kappa0");
		if (attrkappa == nullptr) {
			cout << "missing kappa in bendingspring node!" << endl;
		}
		float kappa = stof(attrkappa->value());

		rapidxml::xml_attribute<>* attrtheta = n->first_attribute("theta0");
		if (attrtheta == nullptr) {
			cout << "missing theta0 in bendingspring node!" << endl;
		}
		float theta = stof(attrtheta->value());

		Particle* p0 = system->particles[p0id];
		Particle* p1 = system->particles[p1id];
		Particle* p2 = system->particles[p2id];
		BendingSpring* s = new BendingSpring(p0, p1, p2);
		// s->kappa0 = kappa;
		s->theta0 = theta;
		return s;
	}
};