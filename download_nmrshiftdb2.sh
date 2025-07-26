#!/bin/bash

# Script to download nmrshiftdb2.xml from SourceForge
# Usage: ./download_nmrshiftdb_xml.sh

set -e

# Try multiple possible XML file URLs
declare -a XML_URLS=(
    "https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2_3d.xml.gz/download"
    "https://sourceforge.net/projects/nmrshiftdb2/files/data/nmrshiftdb2_3d.xml/download"
)

OUTPUT_FILE=""
TEMP_FILE=""

echo "Searching for nmrshiftdb XML files..."
echo ""

# Function to attempt download
attempt_download() {
    local url=$1
    local filename=$(basename $(echo $url | sed 's|/download||'))
    local temp_file="${filename}.tmp"
    
    echo "Trying: $url"
    echo "Output: $filename"
    
    if command -v wget >/dev/null 2>&1; then
        if wget --spider "$url" 2>/dev/null; then
            echo "✓ URL accessible, downloading..."
            wget --progress=bar:force:noscroll \
                 --timeout=30 \
                 --tries=2 \
                 -O "$temp_file" \
                 "$url"
            
            if [ -f "$temp_file" ] && [ -s "$temp_file" ]; then
                mv "$temp_file" "$filename"
                OUTPUT_FILE="$filename"
                return 0
            fi
        else
            echo "✗ URL not accessible"
        fi
    elif command -v curl >/dev/null 2>&1; then
        if curl --head --fail --silent "$url" >/dev/null 2>&1; then
            echo "✓ URL accessible, downloading..."
            curl -L \
                 --progress-bar \
                 --retry 2 \
                 --connect-timeout 30 \
                 -o "$temp_file" \
                 "$url"
            
            if [ -f "$temp_file" ] && [ -s "$temp_file" ]; then
                mv "$temp_file" "$filename"
                OUTPUT_FILE="$filename"
                return 0
            fi
        else
            echo "✗ URL not accessible"
        fi
    fi
    
    # Clean up temp file if it exists
    [ -f "$temp_file" ] && rm "$temp_file"
    return 1
}

# Try each URL
for url in "${XML_URLS[@]}"; do
    if attempt_download "$url"; then
        break
    fi
    echo ""
done

if [ -z "$OUTPUT_FILE" ]; then
    echo "Failed to download XML file from any URL."
    echo ""
    echo "Manual alternatives:"
    echo "1. Visit: https://sourceforge.net/projects/nmrshiftdb2/files/data/"
    echo "2. Look for XML files and download manually"
    echo "3. Or try the SQL approach instead"
    exit 1
fi

echo ""
echo "Download successful!"

# Show file info
FILE_SIZE=$(ls -lh "$OUTPUT_FILE" | awk '{print $5}')
echo "File: $OUTPUT_FILE"
echo "Size: $FILE_SIZE"

# Check if it's compressed
if [[ "$OUTPUT_FILE" == *.gz ]]; then
    echo "File is gzip compressed."
    read -p "Extract now? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Extracting..."
        gunzip "$OUTPUT_FILE"
        EXTRACTED_FILE="${OUTPUT_FILE%.gz}"
        EXTRACTED_SIZE=$(ls -lh "$EXTRACTED_FILE" | awk '{print $5}')
        echo "Extracted: $EXTRACTED_FILE (Size: $EXTRACTED_SIZE)"
        OUTPUT_FILE="$EXTRACTED_FILE"
    fi
fi

# Validate XML and show structure
if [[ "$OUTPUT_FILE" == *.xml ]]; then
    echo ""
    echo "Validating XML structure..."
    
    # Check if it's valid XML
    if python3 -c "import xml.etree.ElementTree as ET; ET.parse('$OUTPUT_FILE')" 2>/dev/null; then
        echo "✓ Valid XML file"
        
        # Show basic structure
        echo ""
        echo "XML structure preview:"
        head -20 "$OUTPUT_FILE" | grep -E "<[^/!?].*?>" | head -10
        
        echo ""
        echo "Root element and first few child elements:"
        python3 -c "
import xml.etree.ElementTree as ET
tree = ET.parse('$OUTPUT_FILE')
root = tree.getroot()
print(f'Root: <{root.tag}>')
for i, child in enumerate(root):
    if i < 5:
        print(f'  Child {i+1}: <{child.tag}>')
    elif i == 5:
        print(f'  ... and {len(root)-5} more children')
        break
"
        
    else
        echo "WARN: File may not be valid XML or may be very large"
    fi
fi

echo ""
echo "Next step:"
echo "Run create_nmrshiftdb2_data.py."
echo ""