#include "InputData.h"

using namespace std;

InputData::InputData():namelist("") {
}


bool BothAreSpaces(char lhs, char rhs) { 
	return (lhs == rhs) && (lhs == ' '); 
}

string InputData::cleanString(string str) {
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	const string whiteSpaces( " \f\n\r\t\v" );
	size_t pos = str.find_last_not_of( whiteSpaces );
	str.erase( pos + 1 );    
	pos = str.find_first_not_of( whiteSpaces );
	str.erase( 0, pos );
	
	std::string::iterator new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
	str.erase(new_end, str.end());   

	return str;
}

void InputData::write(string filename) {
	ofstream ostr(filename.c_str());
	write(ostr);
}

void InputData::write(ostream &ostr) {
	if (ostr!=cerr) ostr << "# smilei " << __VERSION << endl << endl;
	vector<pair<string , map <string,string> > >::iterator  it_type;
	for(it_type = map_.begin(); it_type != map_.end(); it_type++) {
		if (!it_type->first.empty()) ostr << it_type->first << endl;
		map<string, string>::iterator it_type2;
		for(it_type2 = it_type->second.begin(); it_type2 != it_type->second.end(); it_type2++) {
			if (!it_type->first.empty()) ostr << "\t";
			ostr << it_type2->first << " = " << it_type2->second << endl;
		}
		if (!it_type->first.empty()) ostr << "end" << endl;
		ostr << endl;
	}
}

void InputData::parseStream() {
	if (namelist.empty()) ERROR("namelist is empty");
	
	stringstream my_stream(namelist);
	map_.clear();
	map <string,string> defaultmap;
	map <string,string> groupmap;
	
	string strLine ="";
	string group("");

	while (getline(my_stream, strLine)) {
		strLine=strLine.substr(0, strLine.find('#'));				
		strLine=cleanString(strLine);
		if (!strLine.empty()) {
			if (strLine.find('=') == string::npos) {
				if (strLine == "end") {
					map_.push_back(make_pair(group,groupmap));
					group="";
					groupmap.clear();
				} else {
					group=strLine;
				}
			} else {
				stringstream ss(strLine);
				string item;
				while(getline(ss, item, ',')) {
					item=cleanString(item);
					size_t posEqual=item.find('=');
					string left=cleanString(item.substr(0,posEqual));
					string right=cleanString(item.substr(posEqual+1));
					if (group.empty()) {
						defaultmap[left]=right;
					} else {
						groupmap[left]=right;
					}
				}	
			}
		}
	}
	if (!group.empty()) ERROR("Final group "<< group << " not closed. Check the namelist");
	map_.push_back(make_pair("",defaultmap));
}


void InputData::parseFile(string filename) {
	
	ifstream istr(filename.c_str());
	map_.clear();
	map <string,string> defaultmap;
	map <string,string> groupmap;
	
	string strLine ="";
	namelist.clear();
	
	if (istr.is_open()) {
		while (getline(istr, strLine)) {
			strLine=strLine.substr(0, strLine.find('#'));				
			strLine=cleanString(strLine);
			if (!strLine.empty()) namelist += strLine + "\n";
		}
	} else {
		ERROR("File " << filename << "does not exists");
	}
	
}

vector<string> InputData::getGroups() {
	vector<string> vecOut;
	for(vector<pair<string , map <string,string> > >::iterator  it_type = map_.begin(); it_type != map_.end(); it_type++) {
		vecOut.push_back(it_type->first);
	}
	return vecOut;
}

bool InputData::existKey(string key, string group, unsigned int occurrence) {
	unsigned int n_occur=0;
	for (vector<pair<string , map <string,string> > >::iterator  it_type = map_.begin(); it_type != map_.end(); it_type++) {
		if (group == it_type->first) {
			if (occurrence==n_occur) {
				for (map <string,string>::iterator  it_type2 = it_type->second.begin(); it_type2 != it_type->second.end(); it_type2++) {
					if (it_type2->first==key) {
						DEBUG("FOUND KEY");
						return true;
					}
				}
			}
			n_occur++;
		}		
	}
	return true;
}

bool InputData::existGroup(string group, unsigned int occurrence) {
	unsigned int n_occur=0;
	for (vector<pair<string , map <string,string> > >::iterator  it_type = map_.begin(); it_type != map_.end(); it_type++) {
		if (group == it_type->first) {
			if (occurrence==n_occur) {
				return true;
			}
			n_occur++;
		}		
	}
	return false;
}

