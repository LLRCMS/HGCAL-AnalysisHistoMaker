/**
   Editted By     On
   Michael Case   Sun Nov 13 2005
 **/
#include "AnHiMaHGCAL/HGCALGeometry/interface/GeometryConfigurationHack.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DetectorDescription/Parser/interface/DDLParser.h"
#include "DetectorDescription/Base/interface/DDdebug.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>
#include <vector>

GeometryConfigurationHack::GeometryConfigurationHack( const edm::ParameterSet& pset ) : dummyLocation_("") { 
  relFiles_ = pset.getParameter<std::vector<std::string> >("geomXMLFiles");
  for (std::vector<std::string>::const_iterator rit = relFiles_.begin(), ritEnd = relFiles_.end();
      rit != ritEnd; ++rit ) {
    edm::FileInPath fp(*rit);
    files_.push_back(fp.fullPath());
    emptyStrings_.push_back("");
  }
}

GeometryConfigurationHack::~GeometryConfigurationHack() { }

/// Return the Schema Location.
std::string GeometryConfigurationHack::getSchemaLocation() const {
  edm::LogError("GeometryConfigurationHack") << " This sub-class of DDLDocumentProvider does not USE XML parsing!!!" << std::endl;
  return dummyLocation_;
}

/// Return a flag whether to do xml validation or not.
bool GeometryConfigurationHack::doValidation() const {
  LogDebug("GeometryConfigurationHack") << " the doValidation() method not valid for this DDLDocumentProvider" << std::endl;
  return false;
}

/// Return a list of files as a vector of strings.
const std::vector < std::string >  & GeometryConfigurationHack::getFileList(void) const {
  return files_;
}

/// Return a list of urls as a vector of strings.
/**
   The EDM should not allow URLs because of provenance.
   This vector will always be empty.
**/
const std::vector < std::string >  & GeometryConfigurationHack::getURLList(void) const
{
  LogDebug("GeometryConfigurationHack") << " the getURLList of this DDLDocumentProvider empty strings" << std::endl;
  //  return relFiles_;
  return emptyStrings_;
}

/// Print out the list of files.
void GeometryConfigurationHack::dumpFileList(void) const {
  std::cout << "File List:" << std::endl;
  std::cout << "  number of files=" << files_.size() << std::endl;
  for (std::vector<std::string>::const_iterator it = files_.begin(), itEnd = files_.end(); it != itEnd; ++it)
    std::cout << *it << std::endl;
}

int GeometryConfigurationHack::readConfig( const std::string& fname ) {
  edm::LogWarning("GeometryConfigurationHack") << " The readConfig of this DDLDocumentProvider is not valid!" << std::endl;
  return 0;
}

