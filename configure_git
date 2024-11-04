import subprocess

def configure_git():
    email = "hpajoseba@gmail.com"
    username = "jose-bar"
    
    try:
        # Set the email
        subprocess.run(["git", "config", "--global", "user.email", email], check=True)
        # Set the username
        subprocess.run(["git", "config", "--global", "user.name", username], check=True)
        print("Git configuration updated successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")

# Run the function
configure_git()

