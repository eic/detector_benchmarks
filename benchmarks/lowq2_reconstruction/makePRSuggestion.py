#!/usr/bin/env python3
"""
Script to create GitHub PR suggestions for ONNX model updates
"""
import argparse
import sys
from PRfunctions import *

# Parse arguments
parser = argparse.ArgumentParser(description='Update the calibration link for a PR')
parser.add_argument('--pr', type=str, required=True, help='Pull request number')
parser.add_argument('--newURL', type=str, required=True, help='URL of the new updated calibration')
parser.add_argument('--githubToken', type=str, required=True, help='GitHub token for authentication')
parser.add_argument('--calibrationFile', type=str, default='calibrations/onnx/Low-Q2_Steering_Reconstruction.onnx', help='Path to the local calibration file')
parser.add_argument('--xml', type=str, default='compact/calibrations.xml', help='Path to the XML configuration file')
parser.add_argument('--repository', type=str, default='eic/epic', help='GitHub repository (owner/name)')

args = parser.parse_args()


pr_number            = args.pr
new_url              = args.newURL
github_token         = args.githubToken
calibration_file     = args.calibrationFile
xml_file             = args.xml
repository           = args.repository

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

# Create the PR review with suggestion (includes prep if needed)
response = create_pr_suggestion_with_prep(repo_owner, repo_name, pr_number, calibration_file, xml_file, new_url, github_token)

if response:
    print("🎉 PR suggestion completed successfully!")
else:
    print("❌ Failed to create PR suggestion")
    sys.exit(1)

# =============================================================================