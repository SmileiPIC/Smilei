#include "InputData.h"

using namespace std;

string InputData::cleanString(string str) {
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	const string whiteSpaces( " \f\n\r\t\v" );
	size_t pos = str.find_last_not_of( whiteSpaces );
	str.erase( pos + 1 );    
	pos = str.find_first_not_of( whiteSpaces );
	str.erase( 0, pos );
	return str;
}

void InputData::write(string filename, string comment) {
	ofstream ostr(filename.c_str());
	write(ostr, comment);
}

void InputData::write(ostream &ostr, string comment) {
	if (!comment.empty()) ostr << "# " << comment << endl << endl;
	map<string, map<string, string> >::iterator it_type;
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

void InputData::parseFile(string filename) {
	ifstream istr(filename.c_str());
	map_.clear();
	if (istr.is_open()) {
		string strLine ="";
		string group("");
		while (getline(istr, strLine)) {
			strLine=strLine.substr(0, strLine.find('#'));				
			strLine=cleanString(strLine);
			if (!strLine.empty()) {
				if (strLine.find('=') == string::npos) {
					group=strLine;
					if (group == "end") group ="";
				} else {
                    stringstream ss(strLine);
                    string item;
                    while(getline(ss, item, ',')) {
                        item=cleanString(item);
                        size_t posEqual=item.find('=');
                        string left=cleanString(item.substr(0,posEqual));
                        string right=cleanString(item.substr(posEqual+1));
                        map_[group][left]=right;
                    }	
				}
			}
		}
		istr.close();
	} else {
		ERROR("File " << filename << "does not exists");
	}
}

vector<string> InputData::getGroups() {
	vector<string> vecOut;
	map<string, map<string, string> >::iterator it_type;
	for(it_type = map_.begin(); it_type != map_.end(); it_type++) {
		vecOut.push_back(it_type->first);
	}
	return vecOut;
}


