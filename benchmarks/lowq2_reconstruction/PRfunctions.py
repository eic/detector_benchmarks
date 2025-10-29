#!/usr/bin/env python3
"""
Script to create GitHub PR suggestions for ONNX model updates
"""
import requests
import base64
import xml.etree.ElementTree as ET
import os
import mimetypes
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

def upload_image_to_github(repo_owner, repo_name, image_path, github_token):
    """
    Process an image file for use in GitHub comments.
    
    Since GitHub doesn't provide a direct API for uploading images to comments,
    we embed them as base64-encoded data URIs directly in the markdown.
    This works well for reasonably-sized images (< 500KB).
    
    For larger images, they should be uploaded to a separate hosting service
    or committed to the repository and referenced via raw GitHub URLs.
    """
    if not os.path.exists(image_path):
        print(f"❌ Image file not found: {image_path}")
        return None
    
    print(f"Processing image: {image_path}...")
    
    # Read the image file
    try:
        with open(image_path, 'rb') as f:
            image_data = f.read()
    except Exception as e:
        print(f"❌ Failed to read image file: {e}")
        return None
    
    # Detect MIME type
    mime_type, _ = mimetypes.guess_type(image_path)
    if not mime_type or not mime_type.startswith('image/'):
        mime_type = 'application/octet-stream'
    
    # Use data URI to embed the image directly in markdown
    return create_data_uri(image_path, image_data, mime_type)

def create_data_uri(image_path, image_data, mime_type):
    """Create a data URI for embedding images directly in markdown"""
    # Check file size (GitHub has limits on comment size)
    size_kb = len(image_data) / 1024
    if size_kb > 500:  # Limit to 500KB for data URIs
        print(f"⚠️ Image too large for data URI ({size_kb:.1f}KB): {image_path}")
        return None
    
    encoded = base64.b64encode(image_data).decode('utf-8')
    data_uri = f"data:{mime_type};base64,{encoded}"
    print(f"✅ Created data URI for {image_path} ({size_kb:.1f}KB)")
    return data_uri

def process_image_list(image_list, repo_owner, repo_name, github_token):
    """Process a list of images (URLs or local files) and return valid URLs"""
    if not image_list:
        return []
    
    processed_images = []
    for img in image_list:
        # If it's already a URL, use it as-is
        if img.startswith(('http://', 'https://')):
            processed_images.append(img)
        # If it's a data URI, use it as-is
        elif img.startswith('data:'):
            processed_images.append(img)
        # Otherwise, try to upload it
        else:
            uploaded_url = upload_image_to_github(repo_owner, repo_name, img, github_token)
            if uploaded_url:
                processed_images.append(uploaded_url)
            else:
                print(f"⚠️ Skipping image that couldn't be uploaded: {img}")
    
    return processed_images

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

    # Process image lists - upload local files and get URLs
    print("\nProcessing before images...")
    processed_before_images = process_image_list(before_images, repo_owner, repo_name, github_token)
    
    print("\nProcessing after images...")
    processed_after_images = process_image_list(after_images, repo_owner, repo_name, github_token)

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
    if processed_before_images and len(processed_before_images) > 0:
        comment_body += "\n\n---\n\n### 📊 Before Calibration Update\n\n"
        for i, img_url in enumerate(processed_before_images, 1):
            comment_body += f"![Before Image {i}]({img_url})\n\n"
    
    # Add after images section if provided
    if processed_after_images and len(processed_after_images) > 0:
        comment_body += "\n\n---\n\n### 📈 After Calibration Update\n\n"
        for i, img_url in enumerate(processed_after_images, 1):
            comment_body += f"![After Image {i}]({img_url})\n\n"
    
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
            print(f"❌ Failed to update existing comment: {response.status_code}")
            print(f"   Response: {response.text}")
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
            print(f"❌ Failed to create PR comment: {response.status_code}")
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

def find_existing_bot_comment_general(repo_owner, repo_name, pr_number, bot_comment_base, github_token):
    """Find existing bot comment (general PR comment, not line-specific)"""
    print(f"Checking for existing bot comments in PR #{pr_number}...")
    
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
    
    comments = response.json()
    
    # Look for existing bot comment
    for comment in comments:
        # Check if it's from the bot (contains the bot identifier)
        if bot_comment_base in comment.get('body', ''):
            print(f"✅ Found existing bot comment: {comment['id']}")
            return comment['id']
    
    print("No existing bot comment found")
    return None