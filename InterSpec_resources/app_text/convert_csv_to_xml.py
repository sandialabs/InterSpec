import csv
from xml.etree.ElementTree import Element, SubElement, Comment, tostring, fromstring, ParseError
from xml.dom import minidom
import os

# Define the path to your CSV file
csv_file_path = 'InterSpec_strings.csv'

# Define a directory to store the output XML files
output_dir = 'converted_xml_output'
os.makedirs(output_dir, exist_ok=True)

# A dictionary to hold XML data grouped by file_name
xml_data = {}

# Function to safely parse XHTML content
def parse_xhtml(xhtml_str):
    try:
        # Wrap content in a div to ensure a single root element
        wrapped_content = f"<div>{xhtml_str}</div>"
        return fromstring(wrapped_content)
     
        dom = fromstring(wrapped_content)
        children_txt = ''.join(tostring(child, encoding='unicode') for child in dom)
        if dom.text is not None:
            if "Tool to fit for" in xhtml_str:
                print( f"xhtml_str={xhtml_str}")
                print( f"dom.txt={dom.text}, and child join: {children_txt}")
            return dom.text + children_txt
        else:
            if "Tool to fit for" in xhtml_str:
                print( f"No DOM txt")
            return children_txt
    except ParseError:
        # Handle or log parsing errors here
        return None
    
# Read the CSV file
with open(csv_file_path, mode='r', newline='', encoding='utf-8-sig') as file:
    reader = csv.DictReader(file)
    
    # Process each row in the CSV
    for row in reader:
        file_name = row['XML Filename']
        if file_name not in xml_data:
            xml_data[file_name] = []
        xml_data[file_name].append(row)

# Generate XML files
for file_name, messages in xml_data.items():
    # Create the root element
    root = Element('messages')
    
    for msg in messages:
        # Create a message element
        if msg['id']:
            message_elem = SubElement(root, 'message')
            message_elem.set('id', msg['id'])
            if msg['Comment']:
                message_elem.set('comment', msg['Comment'])
            # Parse and append XHTML content
            xhtml_content = parse_xhtml(msg['English Value'])
            #print( tostring(xhtml_content, 'utf-8').decode("utf-8") )
            if xhtml_content is not None:
                message_elem.text = xhtml_content.text
                # Since we wrapped content in a div, we iterate over its children
                for child in xhtml_content:
                    message_elem.append(child)
            #message_elem.text = xhtml_content
            #message_elem.text = msg['English Value']
        else:
            # Add XML comment if present
            #if msg['Comment']:
            #    root.append(Comment(msg['English Value']))
            # comments dont seem to be coming through very well, so we'll do a hack
            message_elem = SubElement(root, 'CommentPlaceholder')
            message_elem.text = msg['English Value']
    
    # Convert the ElementTree to a string
    xml_str = tostring(root, 'utf-8', "xml")
    pretty_xml_str = minidom.parseString(xml_str).toprettyxml(indent="  ")
    
    pretty_xml_str = pretty_xml_str.replace('<CommentPlaceholder>', '<!--' )
    pretty_xml_str = pretty_xml_str.replace('</CommentPlaceholder>', '-->' )
    pretty_xml_str = pretty_xml_str.replace('&quot;', '"' )


    # Write to an XML file
    with open(os.path.join(output_dir, file_name), 'w', encoding='utf-8') as xml_file:
        xml_file.write(pretty_xml_str)
        #xml_file.write(xml_str.decode("utf-8"))
        

print("CSV has been converted back into XML documents.")