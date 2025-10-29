#!/usr/bin/env python3
"""
Script to create GitHub PR comments for ONNX model updates
"""
import argparse
import sys
from PRfunctions import *

# Parse arguments
parser = argparse.ArgumentParser(description='Create a PR comment with calibration update suggestions')
parser.add_argument('--pr', type=str, required=True, help='Pull request number')
parser.add_argument('--newURL', type=str, required=True, help='URL of the new updated calibration')
parser.add_argument('--githubToken', type=str, required=True, help='GitHub token for authentication')
parser.add_argument('--calibrationFile', type=str, default='calibrations/onnx/Low-Q2_Steering_Reconstruction.onnx', help='Path to the local calibration file')
parser.add_argument('--xml', type=str, default='compact/calibrations.xml', help='Path to the XML configuration file')
parser.add_argument('--repository', type=str, default='eic/epic', help='GitHub repository (owner/name)')
parser.add_argument('--beforeImages', type=str, nargs='*', default=[], help='List of before image URLs or paths')
parser.add_argument('--afterImages', type=str, nargs='*', default=[], help='List of after image URLs or paths')


args = parser.parse_args()

pr_number            = args.pr
new_url              = args.newURL
github_token         = args.githubToken
calibration_file     = args.calibrationFile
xml_file             = args.xml
repository           = args.repository
before_images        = args.beforeImages
after_images         = args.afterImages

# =============================================================================

# Extract repo owner and name
repo_owner, repo_name = parse_repository(repository)
if not repo_owner or not repo_name:
    print("Invalid repository format. Exiting.")
    sys.exit(1)


# =============================================================================

print(f"Starting PR update for PR #{pr_number} in {repository}")
pr_info = get_pr_info(repo_owner, repo_name, pr_number, github_token)
if not pr_info:
    print("Failed to retrieve PR information. Exiting.")
    sys.exit(1)

# =============================================================================

content = get_file_content(repo_owner, repo_name, xml_file, pr_info['head']['sha'], github_token)
if not content:
    print("Failed to retrieve file content. Exiting.")
    sys.exit(1)

# =============================================================================

# Parse the XML content and find the line to update
line_number, suggested_line = find_and_update_epic_fileloader_url(content, calibration_file, new_url)

# =============================================================================

if line_number is not None and suggested_line is not None:
    print(f"✅ Found URL to update in {xml_file} at line {line_number}")
    print(f"   Suggested change: {suggested_line.strip()}")
    
    # Create the PR comment with proposed changes
    response = create_pr_suggestion(repo_owner, repo_name, pr_number, calibration_file, xml_file, line_number, suggested_line, pr_info['head']['sha'], github_token, before_images, after_images)
    
    if response:
        print("🎉 PR comment created successfully!")
    else:
        print("❌ Failed to create PR comment")
        sys.exit(1)
else:
    print(f"❌ Failed to find URL to update in {xml_file}")
    sys.exit(1)

# =============================================================================