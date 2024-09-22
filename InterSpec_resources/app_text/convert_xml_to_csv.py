import csv
import xml.etree.ElementTree as ET
import os
import glob

# Define the path for the output CSV file
csv_file_path = 'InterSpec_strings.csv'

class CommentedTreeBuilder(ET.TreeBuilder):
    def comment(self, data):
        self.start(ET.Comment, {})
        self.data(data)
        self.end(ET.Comment)

# Function to extract XHTML content from an element
def get_xhtml_content(element):
    # Serialize the element's children to a string and decode bytes to string
    return ''.join(ET.tostring(child, encoding='unicode') for child in element)


# Open a file for writing with UTF-8 encoding and include the BOM
with open(csv_file_path, mode='w', newline='', encoding='utf-8-sig') as file:
    writer = csv.writer(file)
    
    # Write the header row
    writer.writerow(['XML Filename', 'id', 'Comment', 'English Value'])

    # We could loop over files in a whatever order glob.glob('*.xml') gives us
    #  But lets instead put the files in an order so there is a tiny bit more context
    xml_files = [ 'InterSpecApp.xml', 'InterSpec.xml', 'D3SpectrumDisplayDiv.xml', 'CompactFileManager.xml', \
                 'ReferencePhotopeakDisplay.xml', 'PeakInfoDisplay.xml', 'PeakModel.xml', 'EnergyCalTool.xml', \
                'IsotopeSearchByEnergy.xml', 'DrfSelect.xml', 'ShieldingSelect.xml', \
                'ShieldingSourceDisplay.xml', 'GammaCountDialog.xml', 'MultimediaDisplay.xml', \
                'MakeFwhmForDrf.xml', 'SearchMode3DChart.xml', 'SpecMeasManager.xml', \
                'FeatureMarkerWidget.xml', 'SpecFileSummary.xml', 'LeafletRadMap.xml', \
                'DecayActivity.xml', 'GammaXsGui.xml', 'HelpSystem.xml', 'OneOverR2Calc.xml', 'RemoteRid.xml', \
                 'MakeDrf.xml', 'UnitsConverterTool.xml', 'FluxTool.xml', 'DoseCalcWidget.xml', \
                'ShowRiidInstrumentsAna.xml', 'RelActManualGui.xml', 'LicenseAndDisclaimersWindow.xml' \
                ]
    
    # But lets check to make sure the above array contains all XML files we want

    # Loop over all XML files in the current directory
    for xml_file_path in glob.glob('*.xml'):
        # Extract the file name without the path
        xml_file_name = os.path.basename(xml_file_path)

        # Skip files with "_" in their names
        if "_" in xml_file_name:
            continue

        if xml_file_name not in xml_files:
            print( "Missing file ", xml_file_name, " from `xml_files` array"  )
            exit(-1)

    # Now we can start processing XML file.
    for xml_file_path in xml_files:
        xml_file_name = xml_file_path
        print( "Working on ", xml_file_name )

        # Parse the XML file
        try:
            parser = ET.XMLParser(target=CommentedTreeBuilder())
            tree = ET.parse(xml_file_path, parser)
        except ET.ParseError as e:
            print(f"Error parsing {xml_file_name}: {e}")
            continue  # Skip to the next file if an error occurs
        
        root = tree.getroot()

        # Iterate over each node in the XML
        for node in root:
            # Check if the node is a 'message' element
            if node.tag == 'message':
                # Extract the 'id' attribute, text content of the message, and 'comment' attribute if it exists
                id_attr = node.get('id')
                
                message_text = node.text
                
                # Extract XHTML content
                message_xhtml = get_xhtml_content(node)
                if len(message_xhtml) > 1:
                    if not message_text:
                        message_text = ''
                    message_text = message_text + message_xhtml

                element_comment = node.get('comment', '')  # Default to empty string if 'comment' attribute does not exist
                
                # Write the extracted data to the CSV file, including the element's 'comment' attribute and the last seen XML comment
                writer.writerow([xml_file_name, id_attr, element_comment, message_text])
            
            # Check if the node is a comment
            if node.tag is ET.Comment:
                # Write the comment in the next line with empty values for other columns, except for the comment column
                writer.writerow([xml_file_name, '', '', node.text])

print("Conversion to CSV completed successfully.")
