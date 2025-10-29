#!/usr/bin/env python3
"""
Script to create GitHub PR suggestions for ONNX model updates
"""
import requests
import base64
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime, timezone
import os

# =============================================================================
# Utility Functions
# =============================================================================

def parse_repository(repository):
    """Parse repository string into owner and name"""
    try:
        owner, name = repository.split('/')
        return owner, name
    except ValueError:
        print(f"❌ Invalid repository format: {repository}. Expected: 'owner/name'")
        return None, None

# =============================================================================
# GitHub API Functions
# =============================================================================

def get_pr_info(repo_owner, repo_name, pr_number, github_token=None):
    """Get PR information including head SHA"""
    print(f"Getting PR information for #{pr_number}...")
    
    headers = {
        'Accept': 'application/vnd.github+json'
    }
    
    if github_token:
        headers['Authorization'] = f'token {github_token}'
    
    url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/pulls/{pr_number}"
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        pr_data = response.json()
        print(f"✅ Found PR: {pr_data['title']}")
        print(f"   Branch: {pr_data['head']['ref']}")
        print(f"   SHA: {pr_data['head']['sha']}")
        return pr_data
    else:
        print(f"❌ Failed to get PR info: {response.status_code}")
        print(f"   Response: {response.text}")
        return None

def get_file_content(repo_owner, repo_name, file_path, sha, github_token=None):
    """Get file content from GitHub at a specific commit"""
    print(f"Getting content for {file_path} at commit {sha[:8]}...")
    
    headers = {
        'Accept': 'application/vnd.github+json'
    }
    
    if github_token:
        headers['Authorization'] = f'token {github_token}'
    
    url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{file_path}?ref={sha}"
    response = requests.get(url, headers=headers)
    
    if response.status_code == 200:
        content_data = response.json()
        if 'content' in content_data:
            # GitHub returns base64-encoded content
            content = base64.b64decode(content_data['content']).decode('utf-8')
            print(f"✅ Got file content ({len(content)} characters)")
            return content
        else:
            print(f"❌ No content found in response")
            return None
    elif response.status_code == 404:
        print(f"❌ File not found: {file_path}")
        return None
    else:
        print(f"❌ Failed to get file content: {response.status_code}")
        print(f"   Response: {response.text}")
        return None
    
# =============================================================================
# XML Processing Functions
# =============================================================================

def find_and_update_epic_fileloader_url(content, calibration_file, new_url, plugin_name='epic_FileLoader'):
    """Find the line with FileLoader URL and return line number and suggested change"""
    print(f"Parsing XML to find {plugin_name} URL for {calibration_file}...")
    
    try:
        # Parse the XML content
        root = ET.fromstring(content)
        calibration_filename = Path(calibration_file).name
        
        # Find all plugin elements with the specified name
        file_loader_plugins = root.findall(f".//plugin[@name='{plugin_name}']")
        
        if not file_loader_plugins:
            print(f"❌ No {plugin_name} plugin found")
            return None, None
        
        print(f"✅ Found {len(file_loader_plugins)} {plugin_name} plugin(s)")
        
        # Look through each FileLoader plugin
        for plugin in file_loader_plugins:
            print(f"   Checking plugin...")
            
            # Find all arg elements with value starting with "url:"
            for arg in plugin.findall("arg"):
                value = arg.get('value', '')
                
                if value.startswith('url:') and calibration_filename in value:
                    print(f"✅ Found matching URL arg: {value}")
                    
                    # Find the line number of this change in the original content
                    line_number = find_line_number_of_change(content, value)
                    
                    if line_number:
                        # Get the original line and create the suggested replacement
                        lines = content.split('\n')
                        original_line = lines[line_number - 1]  # Convert to 0-based index
                        suggested_line = original_line.replace(value, f'url:{new_url}')
                        
                        print(f"✅ Line {line_number}: {original_line.strip()}")
                        print(f"   Suggested: {suggested_line.strip()}")
                        
                        return line_number, suggested_line
                    else:
                        print("❌ Could not find line number for the URL")
                        return None, None
        
        print(f"❌ No matching URL argument found in {plugin_name} plugins")
        return None, None
        
    except ET.ParseError as e:
        print(f"❌ XML parsing error: {e}")
        return None, None
    
def find_line_number_of_change(original_content, old_value):
    """Find the line number where the change occurred"""
    lines = original_content.split('\n')
    
    for i, line in enumerate(lines, 1):
        if old_value in line:
            return i
    
    return None

# =============================================================================
# GitHub PR Comment Functions
# =============================================================================

def create_pr_suggestion(repo_owner, repo_name, pr_number, calibration_file, xml_file, line_number, suggested_line, head_sha, github_token, before_images=None, after_images=None):
    """Create a PR comment with proposed changes"""
    print(f"Creating PR comment with calibration update for #{pr_number}...")
    
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }

    bot_comment_base = f"🤖 **Automated Calibration `{calibration_file}` Update**"

    # Check for existing comments from bot
    existing_comment_id = find_existing_bot_comment_general(repo_owner, repo_name, pr_number, bot_comment_base, github_token)
    
    # Generate timestamp for the comment
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    # Get current line content for context
    content = get_file_content(repo_owner, repo_name, xml_file, head_sha, github_token)
    lines = content.split('\n') if content else []
    current_line = lines[line_number - 1].strip() if line_number <= len(lines) else "Line not found"

    comment_body = f"""{bot_comment_base}{' (Updated)' if existing_comment_id else ''}

A new calibration has been generated and is ready for use.

**File:** `{xml_file}`
**Line:** {line_number}
**Last updated:** {timestamp}

**Current line:**
```xml
{current_line}
```

**Proposed change:**
```xml
{suggested_line.strip()}
```

Please update the calibration URL in `{xml_file}` at line {line_number}."""

    # Add before images section if provided
    if before_images:
        comment_body += "\n\n---\n\n### 📊 Before Calibration Update\n\n"
        for img_url in before_images:
            comment_body += f"![Before Image]({img_url})\n\n"

    # Add after images section if provided
    if after_images:
        comment_body += "\n\n---\n\n### 📈 After Calibration Update\n\n"
        for img_url in after_images:
            comment_body += f"![After Image]({img_url})\n\n"
    
    if existing_comment_id:
        # Update existing comment
        print(f"Updating existing comment {existing_comment_id}...")
        update_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/issues/comments/{existing_comment_id}"
        update_data = {'body': comment_body}
        response = requests.patch(update_url, headers=headers, json=update_data)
        if response.status_code == 200:
            print("✅ Existing PR comment updated successfully")
            return response.json()
        else:
            print(f"❌ Failed to update existing comment: {response.status_code}\n{response.text}")
            return None
    else:
        # Create new regular PR comment
        print("Creating new PR comment...")
        comment_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/issues/{pr_number}/comments"
        comment_data = {'body': comment_body}
        response = requests.post(comment_url, headers=headers, json=comment_data)
        if response.status_code == 201:
            print("✅ New PR comment created successfully")
            return response.json()
        else:
            print(f"❌ Failed to create PR comment: {response.status_code}\n{response.text}")
            return None

def find_existing_bot_comment_general(repo_owner, repo_name, pr_number, bot_comment_base, github_token):
    """Find existing bot comment (general PR comment, not line-specific)"""
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }
        
    # Get all general comments for the PR
    comments_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/issues/{pr_number}/comments"
    response = requests.get(comments_url, headers=headers)
    
    if response.status_code != 200:
        print(f"❌ Failed to get PR comments: {response.status_code}")
        return None
    
    # Look for existing bot comment
    for comment in response.json():
        if bot_comment_base in comment.get('body', ''):
            print(f"✅ Found existing bot comment: {comment['id']}")
            return comment['id']
    
    print("No existing bot comment found")
    return None