import React from "react";

interface InputProps {
  type?: "file"; // make optional, defaults to file
  value?: File | null; // optional now
  onChange: (event: React.ChangeEvent<HTMLInputElement>) => void;
  placeholder?: string;
  className?: string;
  disabled?: boolean;
  accept?: string;
}

const Input: React.FC<InputProps> = ({
  type = "file",
  value,
  onChange,
  placeholder,
  className = "",
  disabled = false,
  accept,
}) => {
  return (
    <div className={`input-wrapper ${className}`}>
      <input
        type={type}
        // ⚠️ do not bind value for file inputs
        onChange={onChange}
        placeholder={placeholder}
        className="input"
        disabled={disabled}
        accept={accept}
      />
      {/* Display the name of the selected file */}
      {value && <p>Selected file: {value.name}</p>}
    </div>
  );
};

export default Input;
