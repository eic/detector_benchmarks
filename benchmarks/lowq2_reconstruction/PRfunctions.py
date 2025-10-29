#!/usr/bin/env python3
"""
Script to create GitHub PR suggestions for ONNX model updates
"""
import requests
import base64
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime, timezone

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

def make_github_request(url, headers, method='GET', data=None):
    """Make a GitHub API request with consistent error handling"""
    try:
        if method.upper() == 'GET':
            response = requests.get(url, headers=headers)
        elif method.upper() == 'POST':
            response = requests.post(url, headers=headers, json=data)
        elif method.upper() == 'PATCH':
            response = requests.patch(url, headers=headers, json=data)
        else:
            print(f"❌ Unsupported HTTP method: {method}")
            return None
        
        return response
    except requests.RequestException as e:
        print(f"❌ Request failed: {e}")
        return None

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

def create_pr_suggestion(repo_owner, repo_name, pr_number, calibration_file, xml_file, line_number, suggested_line, head_sha, github_token):
    """Create a PR review with code suggestion or update existing one"""
    print(f"Creating/updating PR review with suggestion for #{pr_number}...")
    
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }

    bot_comment_base = f"🤖 **Automated Calibration `{calibration_file}` Update**"

    # Check for existing review comments from bot
    existing_comment_id = find_existing_bot_comment(repo_owner, repo_name, pr_number, bot_comment_base, xml_file, line_number, github_token)
    
    # Generate timestamp for the comment
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    suggestion_body = f"""{bot_comment_base}{' (Updated)' if existing_comment_id else ''}

A new calibration has been generated and is ready for use.

**Last updated:** {timestamp}

```suggestion
{suggested_line}
```"""
    
    if existing_comment_id:
        # Update existing comment
        print(f"Updating existing comment {existing_comment_id}...")
        update_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/pulls/comments/{existing_comment_id}"
        update_data = {'body': suggestion_body}
        response = requests.patch(update_url, headers=headers, json=update_data)
        
        if response.status_code == 200:
            print("✅ Existing PR comment updated successfully")
            return response.json()
        else:
            print(f"❌ Failed to update existing comment: {response.status_code}")
            print(f"   Response: {response.text}")
            return None
    else:
        # Create new review comment
        print("Creating new review comment...")
        review_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/pulls/{pr_number}/reviews"
        
        review_data = {
            'body': f'🤖 Automated review with updated calibration URLs `{calibration_file}` for PR #{pr_number}',
            'event': 'COMMENT',
            'commit_id': head_sha,
            'comments': [
                {
                    'path': xml_file,
                    'line': line_number,
                    'body': suggestion_body
                }
            ]
        }
        
        response = requests.post(review_url, headers=headers, json=review_data)
        
        if response.status_code == 200:
            print("✅ New PR review with suggestion created successfully")
            return response.json()
        else:
            print(f"❌ Failed to create PR review: {response.status_code}")
            print(f"   Response: {response.text}")
            return None

def find_existing_bot_comment(repo_owner, repo_name, pr_number, bot_comment_base, xml_file, line_number, github_token):
    """Find existing bot comment on the specific line"""
    print(f"Checking for existing bot comments on line {line_number}...")
    
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }
        
    # Get all review comments for the PR
    comments_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/pulls/{pr_number}/comments"
    response = requests.get(comments_url, headers=headers)
    
    if response.status_code != 200:
        print(f"❌ Failed to get PR comments: {response.status_code}")
        return None
    
    comments = response.json()
    
    # Look for existing bot comment on the same line and file
    for comment in comments:
        # Check if it's from the bot (contains the bot identifier)
        # and on the same file and line
        if (comment.get('path') == xml_file and 
            comment.get('line') == line_number and
            bot_comment_base in comment.get('body', '')):
            print(f"✅ Found existing bot comment: {comment['id']}")
            return comment['id']
    
    print("No existing bot comment found")
    return None

def create_commit_with_marker_for_suggestion(repo_owner, repo_name, pr_number, xml_file, line_number, github_token):
    """Create a minimal commit that adds the XML file to PR changes to enable suggestions"""
    print(f"Creating minimal commit to enable suggestions on {xml_file}...")
    
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }
    
    # Get PR info to get the branch
    pr_info = get_pr_info(repo_owner, repo_name, pr_number, github_token)
    if not pr_info:
        return False
        
    branch_name = pr_info['head']['ref']
    head_sha = pr_info['head']['sha']
    
    print(f"Target branch: {branch_name}")
    print(f"Current head: {head_sha}")
    
    # Get current file content
    file_content = get_file_content(repo_owner, repo_name, xml_file, head_sha, github_token)
    if not file_content:
        print(f"❌ Could not get content for {xml_file}")
        return False
    
    # Get the current file SHA (required for updates)
    file_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{xml_file}?ref={head_sha}"
    response = requests.get(file_url, headers=headers)
    
    if response.status_code != 200:
        print(f"❌ Could not get file SHA: {response.status_code}")
        return False
        
    file_sha = response.json().get('sha')
    
    # Add a minimal comment to the line that will be suggested
    lines = file_content.split('\n')
    if line_number <= len(lines):
        # Add a comment at the end of the target line
        target_line = lines[line_number - 1]  # Convert to 0-based index
        
        # Check if comment already exists
        if '<!-- calibration-update-marker -->' not in target_line:
            lines[line_number - 1] = target_line + ' <!-- calibration-update-marker -->'
            updated_content = '\n'.join(lines)
            
            print(f"Adding marker comment to line {line_number}")
            
            # Create the commit
            commit_data = {
                'message': '[skip ci] Add calibration update marker\n\nMinimal change to enable calibration suggestions.',
                'content': base64.b64encode(updated_content.encode()).decode(),
                'sha': file_sha,
                'branch': branch_name
            }
            
            commit_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{xml_file}"
            response = requests.put(commit_url, headers=headers, json=commit_data)
            
            if response.status_code in [200, 201]:
                print("✅ Minimal commit created successfully")
                return True
            else:
                print(f"❌ Failed to create minimal commit: {response.status_code}")
                print(f"   Response: {response.text}")
                return False
        else:
            print("ℹ️ Marker comment already exists")
            return True
    else:
        print(f"❌ Line number {line_number} is out of range")
        return False

def create_pr_suggestion_with_prep(repo_owner, repo_name, pr_number, calibration_file, xml_file, new_url, github_token):
    """Create PR suggestion with preparation if needed"""
    print(f"Creating PR suggestion with preparation for #{pr_number}...")
    
    # Get current PR info
    pr_info = get_pr_info(repo_owner, repo_name, pr_number, github_token)
    if not pr_info:
        return False
    
    head_sha = pr_info['head']['sha']
    
    # Get file content and find the line to update
    file_content = get_file_content(repo_owner, repo_name, xml_file, head_sha, github_token)
    if not file_content:
        return False
    
    line_number, suggested_line = find_and_update_epic_fileloader_url(file_content, calibration_file, new_url)
    if not line_number or not suggested_line:
        print(f"❌ Could not find calibration line to update in {xml_file}")
        return False
    
    # Check if XML file is already in the PR changes
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }
    
    files_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/pulls/{pr_number}/files"
    response = requests.get(files_url, headers=headers)
    
    if response.status_code == 200:
        changed_files = response.json()
        file_paths = [f['filename'] for f in changed_files]
        
        if xml_file not in file_paths:
            print(f"📝 {xml_file} is not in PR changes, creating minimal commit to enable suggestions...")
            
            # Create minimal commit to add file to PR
            success = create_commit_with_marker_for_suggestion(repo_owner, repo_name, pr_number, xml_file, line_number, github_token)
            
            if not success:
                print("❌ Failed to create minimal commit, falling back to regular comment")
                return create_pr_comment_fallback(repo_owner, repo_name, pr_number, calibration_file, new_url, github_token)
            
            # Wait a moment for GitHub to process the commit
            import time
            time.sleep(3)
            
            # Get updated PR info with new head SHA
            pr_info = get_pr_info(repo_owner, repo_name, pr_number, github_token)
            if not pr_info:
                return False
                
            head_sha = pr_info['head']['sha']
            
            # Get updated file content (now with the marker comment)
            file_content = get_file_content(repo_owner, repo_name, xml_file, head_sha, github_token)
            if not file_content:
                return False
            
            # Re-find the line (it should be the same, but content is slightly different)
            line_number, suggested_line = find_and_update_epic_fileloader_url(file_content, calibration_file, new_url)
            if not line_number or not suggested_line:
                print(f"❌ Could not find calibration line after minimal commit")
                return False
            
            # Remove the marker comment from the suggested line
            suggested_line = suggested_line.replace(' <!-- calibration-update-marker -->', '')
            
            print(f"✅ XML file now in PR changes, proceeding with suggestion")
        else:
            print(f"✅ {xml_file} already in PR changes")
    else:
        print(f"❌ Could not get PR files: {response.status_code}")
        return False
    
    # Now create the proper code suggestion
    return create_pr_suggestion(repo_owner, repo_name, pr_number, calibration_file,
                              xml_file, line_number, suggested_line, head_sha, github_token)

def create_pr_comment_fallback(repo_owner, repo_name, pr_number, calibration_file, new_url, github_token):
    """Fallback: Create a regular PR comment when suggestions aren't possible"""
    print(f"Creating fallback PR comment...")
    
    headers = {
        'Accept': 'application/vnd.github+json',
        'Authorization': f'token {github_token}'
    }
    
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    
    comment_body = f"""🤖 **Automated Calibration `{calibration_file}` Update**

A new calibration has been generated and is ready for use.

**New URL:** `{new_url}`
**Last updated:** {timestamp}

**Manual Update Required:**
Please update the XML configuration to use the new calibration URL.

```xml
url:{new_url}
```

Look for the existing `{Path(calibration_file).name}` entry and update its URL."""
    
    comment_url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/issues/{pr_number}/comments"
    comment_data = {'body': comment_body}
    
    response = requests.post(comment_url, headers=headers, json=comment_data)
    
    if response.status_code == 201:
        print("✅ Fallback PR comment created successfully")
        return response.json()
    else:
        print(f"❌ Failed to create fallback comment: {response.status_code}")
        return None